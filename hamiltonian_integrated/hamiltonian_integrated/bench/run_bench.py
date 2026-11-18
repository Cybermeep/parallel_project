import argparse
import csv
import gzip
import os
import shlex
import statistics
import subprocess
import sys
import time
import urllib.request
from collections import OrderedDict, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

RESULT_FIELDS = ["instance", "mode", "algorithm", "ranks", "threads", "vertices", "edges", "status", "cost",
                 "expected", "verdict", "seconds", "nodes", "lower_bound", "verified"]

OUT_FIELDS = ["suite", "instance", "mode", "algo", "algorithm", "ranks", "threads", "p", "vertices", "edges",
              "status", "cost", "expected", "verdict", "seconds", "min_seconds", "max_seconds", "wall_seconds",
              "repeats", "nodes", "lower_bound", "verified", "T1", "speedup", "efficiency", "serial_fraction",
              "amdahl_bound", "source", "command"]

TSPLIB_URLS = ["http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/{name}.tsp.gz",
               "https://raw.githubusercontent.com/mastqe/tsplib/master/{name}.tsp"]


def parse_list(text, cast=str):
    return [cast(x.strip()) for x in text.split(",") if x.strip()]


def parse_args():
    p = argparse.ArgumentParser(
        description="Run the Hamiltonian path / TSP solver on public benchmark instances, check every answer "
                    "against the published result, and report S(p) = T(1)/T(p), E(p) = S(p)/p and the "
                    "Amdahl serial fraction.")
    p.add_argument("--bin", default=os.path.join(ROOT, "hampath"))
    p.add_argument("--mpi-bin", default=os.path.join(ROOT, "hampath_mpi"))
    p.add_argument("--launcher", default="mpirun -np {ranks}",
                   help="MPI launch prefix; {ranks} and {threads} are substituted")
    p.add_argument("--manifest", action="append", default=None)
    p.add_argument("--suite", default="smoke", help="comma list of suites, or all")
    p.add_argument("--filter", default="", help="only instances whose name contains this text")
    p.add_argument("--algos", default="pruned",
                   help="pruned,sam-iterative,sam-recursive,hcp-as-tsp (graph modes only)")
    p.add_argument("--threads", default="1,2,4", help="OpenMP threads per rank")
    p.add_argument("--ranks", default="1", help="MPI ranks; values above 1 use --mpi-bin")
    p.add_argument("--repeat", type=int, default=1)
    p.add_argument("--time-limit", type=float, default=60.0)
    p.add_argument("--frontier", type=int, default=16)
    p.add_argument("--out", default=os.path.join(ROOT, "results", "bench.csv"))
    p.add_argument("--summary", default="")
    p.add_argument("--baseline", action="append", default=[],
                   help="earlier result CSV(s) to take T(1) from when this run has no p = 1 configuration")
    p.add_argument("--plots", default="")
    p.add_argument("--fetch-tsplib", action="store_true",
                   help="download the TSPLIB files listed in the manifest into data/tsplib and exit")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def load_manifest(paths, suites, text_filter):
    rows = []
    for path in paths:
        with open(path) as f:
            for row in csv.DictReader(f):
                if suites != ["all"] and row["suite"] not in suites:
                    continue
                if text_filter and text_filter not in row["instance"]:
                    continue
                if not os.path.isabs(row["file"]):
                    row["file"] = os.path.join(ROOT, row["file"])
                rows.append(row)
    return rows


def fetch_tsplib(instances):
    dest = os.path.join(ROOT, "data", "tsplib")
    os.makedirs(dest, exist_ok=True)
    ok, failed = 0, []
    for inst in instances:
        if not inst["file"].endswith(".tsp"):
            continue
        name = os.path.basename(inst["file"])[:-4]
        target = os.path.join(dest, name + ".tsp")
        if os.path.exists(target) and os.path.getsize(target) > 0:
            ok += 1
            continue
        done = False
        for pattern in TSPLIB_URLS:
            url = pattern.format(name=name)
            try:
                data = urllib.request.urlopen(url, timeout=60).read()
                if url.endswith(".gz"):
                    data = gzip.decompress(data)
                with open(target, "wb") as f:
                    f.write(data)
                done = True
                break
            except Exception:
                continue
        if done:
            ok += 1
        else:
            failed.append(name)
    print("TSPLIB files ready: {}; failed: {}".format(ok, ", ".join(failed) if failed else "none"))
    return 0 if not failed else 1


def build_command(args, inst, algo, ranks, threads):
    cmd = []
    if ranks > 1:
        cmd += shlex.split(args.launcher.format(ranks=ranks, threads=threads))
        cmd.append(args.mpi_bin)
    else:
        cmd.append(args.bin)
    cmd += [inst["file"], "--mode", inst["mode"], "--threads", str(threads), "--time-limit", str(args.time_limit),
            "--frontier", str(args.frontier), "--label", inst["instance"], "--quiet"]
    batch = "--batch" in inst.get("args", "")
    if inst["mode"] in ("path", "cycle"):
        cmd += ["--algo", algo]
        if algo.startswith("sam") and inst["mode"] == "path" and not batch:
            cmd.append("--all-starts")
    if inst.get("expect"):
        cmd += ["--expect", inst["expect"]]
    if inst.get("args"):
        cmd += shlex.split(inst["args"])
    return cmd


def run_one(cmd, threads, limit):
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    env.setdefault("OMP_PROC_BIND", "close")
    env.setdefault("OMP_PLACES", "cores")
    start = time.time()
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                              env=env, timeout=limit * 2 + 120)
    except subprocess.TimeoutExpired:
        return None, time.time() - start, "killed after hard timeout"
    except OSError as e:
        return None, time.time() - start, str(e)
    wall = time.time() - start
    lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT,")]
    if not lines:
        return None, wall, (proc.stderr or proc.stdout).strip()[-400:]
    values = lines[-1].split(",", 1)[1].split(",")
    return dict(zip(RESULT_FIELDS, values)), wall, ""


def solved(status):
    if status.startswith("BATCH("):
        return "timeout=0" in status
    return status in ("FOUND", "NONE", "OPTIMAL")


def configurations(args, inst, algos):
    batch = "--batch" in inst.get("args", "")
    combos = []
    for algo in (algos if inst["mode"] in ("path", "cycle") else ["bnb"]):
        for r in parse_list(args.ranks, int):
            for t in parse_list(args.threads, int):
                if algo.startswith("sam") and (r > 1 or t > 1) and not batch:
                    continue
                combos.append((algo, r, t))
    return combos


def load_baselines(paths):
    table = {}
    for path in paths:
        if not os.path.exists(path):
            print("baseline file not found:", path)
            continue
        with open(path) as f:
            for row in csv.DictReader(f):
                if row.get("p") == "1" and solved(row.get("status", "")):
                    table[(row["instance"], row["algo"])] = float(row["seconds"])
    return table


def apply_formulas(rows, baselines):
    groups = defaultdict(list)
    for r in rows:
        groups[(r["suite"], r["instance"], r["algo"])].append(r)
    for (suite, instance, algo), items in groups.items():
        base = [r for r in items if r["p"] == 1]
        if base:
            t1 = float(base[0]["seconds"]) if solved(base[0]["status"]) else None
        else:
            t1 = baselines.get((instance, algo))
        for r in items:
            r["T1"] = t1
            r["speedup"] = r["efficiency"] = r["serial_fraction"] = r["amdahl_bound"] = None
            tp = float(r["seconds"])
            if t1 is None or not solved(r["status"]) or tp <= 0:
                continue
            p = r["p"]
            s = t1 / tp
            r["speedup"] = s
            r["efficiency"] = s / p
            if p > 1 and s > 0:
                f = (1.0 / s - 1.0 / p) / (1.0 - 1.0 / p)
                if 0.0 <= f <= 1.0:
                    r["serial_fraction"] = f
                    if f > 0:
                        r["amdahl_bound"] = 1.0 / (f + (1.0 - f) / p)


def fnum(value, digits=2):
    if value is None or value == "":
        return "-"
    return "{:.{}f}".format(value, digits)


def write_csv(path, rows):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_FIELDS)
        w.writeheader()
        for r in rows:
            out = {k: r.get(k, "") for k in OUT_FIELDS}
            for k in ("T1", "speedup", "efficiency", "serial_fraction", "amdahl_bound"):
                out[k] = "" if out[k] is None else "{:.6f}".format(out[k])
            w.writerow(out)


def write_summary(path, rows, args):
    lines = ["# Benchmark summary", ""]
    lines.append("Every answer is checked against the published result for its public instance "
                 "(Hamiltonicity from the Flinders sets, optimal tour length from TSPLIB).")
    lines.append("")
    lines.append("Formulas (IT 388): S(p) = T(1) / T(p), E(p) = S(p) / p, Amdahl S(p) <= 1 / (f + (1 - f) / p). "
                 "p = MPI ranks x OpenMP threads. The serial fraction f is solved from Amdahl's law using the "
                 "measured S(p): f = (1/S - 1/p) / (1 - 1/p); it is shown only when 1 <= S(p) <= p, because a slowdown or a superlinear search speedup has no Amdahl fraction. In MPI runs rank 0 only hands out work, so p counts "
                 "one core per thread on rank 0 that does no searching.")
    lines.append("")
    lines.append("Repeats per configuration: {}; T(p) is the median. Time limit: {} s.".format(args.repeat,
                                                                                              args.time_limit))
    lines.append("")
    lines.append("## Correctness against published answers")
    lines.append("")
    lines.append("| suite | runs | MATCH | MISMATCH | UNSOLVED / GAP / UNKNOWN | INVALID / ERROR |")
    lines.append("|---|---|---|---|---|---|")
    by_suite = OrderedDict()
    for r in rows:
        by_suite.setdefault(r["suite"], []).append(r)
    for suite, items in by_suite.items():
        c = defaultdict(int)
        for r in items:
            v = r["verdict"]
            c["MATCH" if v == "MATCH" else "MISMATCH" if v == "MISMATCH" else "BAD" if v in ("INVALID", "ERROR")
              else "OTHER"] += 1
        lines.append("| {} | {} | {} | {} | {} | {} |".format(suite, len(items), c["MATCH"], c["MISMATCH"],
                                                             c["OTHER"], c["BAD"]))
    lines.append("")
    for suite, items in by_suite.items():
        lines.append("## Suite `{}`".format(suite))
        lines.append("")
        lines.append("| instance | n | published | algo | ranks x threads | p | status | cost | verdict | T(p) s | "
                     "S(p) | E(p) | f | Amdahl max S(p) |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in items:
            cost = r["cost"] if r["cost"] not in ("", "-1") else "-"
            lines.append("| {} | {} | {} | {} | {} x {} | {} | {} | {} | {} | {:.4f} | {} | {} | {} | {} |".format(
                r["instance"], r["vertices"], r["expected"] or "-", r["algo"], r["ranks"], r["threads"], r["p"],
                r["status"], cost, r["verdict"], float(r["seconds"]), fnum(r["speedup"]), fnum(r["efficiency"]),
                fnum(r["serial_fraction"], 3), fnum(r["amdahl_bound"])))
        lines.append("")
    optimal = defaultdict(set)
    for r in rows:
        if r["status"] == "OPTIMAL":
            optimal[r["instance"]].add(r["cost"])
    disagree = sorted(k for k, v in optimal.items() if len(v) > 1)
    lines.append("## Consistency")
    lines.append("")
    lines.append("All configurations agree on every proven optimum." if not disagree
                 else "Configurations disagree on: " + ", ".join(disagree))
    lines.append("")
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(lines))


def make_plots(rows, folder):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plots")
        return
    os.makedirs(folder, exist_ok=True)
    groups = defaultdict(lambda: defaultdict(list))
    for r in rows:
        if r["speedup"] is not None:
            groups[(r["suite"], r["algo"])][r["instance"]].append((r["p"], r["speedup"], r["efficiency"]))
    for (suite, algo), series in groups.items():
        series = {k: sorted(v) for k, v in series.items() if len(v) > 1}
        if not series:
            continue
        pmax = max(p for pts in series.values() for p, _, _ in pts)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
        axes[0].plot([1, pmax], [1, pmax], linestyle="--", color="#999999", linewidth=1, label="ideal S(p) = p")
        axes[1].axhline(1.0, linestyle="--", color="#999999", linewidth=1)
        for name, pts in sorted(series.items()):
            axes[0].plot([p for p, _, _ in pts], [s for _, s, _ in pts], marker="o", linewidth=1.5, label=name)
            axes[1].plot([p for p, _, _ in pts], [e for _, _, e in pts], marker="o", linewidth=1.5, label=name)
        axes[0].set_xlabel("p (ranks x threads)")
        axes[0].set_ylabel("S(p) = T(1) / T(p)")
        axes[1].set_xlabel("p (ranks x threads)")
        axes[1].set_ylabel("E(p) = S(p) / p")
        for ax in axes:
            try:
                ax.set_xscale("log", base=2)
            except TypeError:
                ax.set_xscale("log", basex=2)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        fig.suptitle("{} - {}".format(suite, algo))
        if len(series) <= 12:
            axes[0].legend(fontsize=7, frameon=False)
        fig.tight_layout()
        path = os.path.join(folder, "scaling_{}_{}.png".format(suite, algo))
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print("plot", path)


def main():
    args = parse_args()
    manifests = args.manifest or [os.path.join(ROOT, "bench", "instances.csv")]
    instances = load_manifest(manifests, parse_list(args.suite), args.filter)
    if args.fetch_tsplib:
        return fetch_tsplib(instances)
    if not instances:
        print("no instances selected")
        return 1
    algos = parse_list(args.algos)
    rows, missing = [], []
    for inst in instances:
        if not os.path.exists(inst["file"]):
            missing.append(inst["instance"])
            continue
        for algo, ranks, threads in configurations(args, inst, algos):
            cmd = build_command(args, inst, algo, ranks, threads)
            if args.dry_run:
                print(" ".join(shlex.quote(c) for c in cmd))
                continue
            runs = []
            for _ in range(args.repeat):
                res, wall, err = run_one(cmd, threads, args.time_limit)
                if res is None:
                    print("  error:", err)
                    res = {k: "" for k in RESULT_FIELDS}
                    res.update(status="ERROR", verdict="ERROR", seconds=str(wall), cost="-1")
                res["wall_seconds"] = "{:.4f}".format(wall)
                runs.append(res)
            times = [float(r["seconds"] or 0) for r in runs]
            mid = sorted(runs, key=lambda r: float(r["seconds"] or 0))[len(runs) // 2]
            row = {k: mid.get(k, "") for k in RESULT_FIELDS}
            row.update(suite=inst["suite"], instance=inst["instance"], mode=inst["mode"], algo=algo, ranks=ranks,
                       threads=threads, p=ranks * threads, expected=inst.get("expect", ""),
                       seconds="{:.6f}".format(statistics.median(times)), min_seconds="{:.6f}".format(min(times)),
                       max_seconds="{:.6f}".format(max(times)), wall_seconds=mid["wall_seconds"],
                       repeats=len(runs), source=inst.get("source", ""),
                       command=" ".join(shlex.quote(c) for c in cmd))
            rows.append(row)
            print("{:10}{:30}{:15}{:>3}x{:<4}{:38}{:11}{:>11.4f}s".format(
                inst["suite"], inst["instance"][:29], algo, ranks, threads, row["status"][:37], row["verdict"],
                float(row["seconds"])))
            sys.stdout.flush()
    if args.dry_run:
        return 0
    if missing:
        print("skipped {} instance(s) with no data file (try --fetch-tsplib on a login node): {}".format(
            len(missing), ", ".join(missing)))
    if not rows:
        return 1
    apply_formulas(rows, load_baselines(args.baseline))
    write_csv(args.out, rows)
    summary = args.summary or os.path.splitext(args.out)[0] + "_summary.md"
    write_summary(summary, rows, args)
    print("wrote", args.out)
    print("wrote", summary)
    if args.plots:
        make_plots(rows, args.plots)
    bad = [r for r in rows if r["verdict"] in ("MISMATCH", "INVALID")]
    if bad:
        print("{} run(s) disagree with the published answer".format(len(bad)))
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
