import re
from pathlib import Path

# --------------------------------------------------
# Configuration
# --------------------------------------------------

LOG_FILE = "benchmark_results.log"
OUTPUT_FILE = "generated_tables.tex"

N_VALUES = [2**13, 2**15, 2**17]
L_VALUES = [8, 16, 24]
THREADS = [1, 2, 4]

# Table 1 has no L dimension, so select one L.
PRIMITIVE_L = 8

# Number of threads used for the OpenFHE comparison.
OPENFHE_COMPARE_THREADS = 4

# --------------------------------------------------
# Regular expressions
# --------------------------------------------------

COMMAND_RE = re.compile(r"Command:\s+make\s+(.*)")
PARAM_RE = re.compile(r"(\w+)=([^\s]+)")

AVG_RE = re.compile(
    r"^(.*?)\)\s+Total:.*?AVG:\s*([0-9.]+)\s+us"
)

# --------------------------------------------------
# Parse benchmark log
# --------------------------------------------------

perf_data = {}
openfhe_data = {}

current = None
successful = True

with open(LOG_FILE, "r") as log:

    for line in log:

        # Identify a new command.
        match = COMMAND_RE.search(line)

        if match:
            params = dict(PARAM_RE.findall(match.group(1)))

            current = {
                "target": "openfhe" if "openfhe" in match.group(1).split()
                           else "perf",
                "N": int(params["N"]),
                "L": int(params["L"]),
                "OPT": int(params["OPT"]) if "OPT" in params else None,
                "OMP": int(params["OMP_NUM_THREADS"]),
                "times": {}
            }

            continue

        if current is None:
            continue

        # Extract execution times.
        match = AVG_RE.search(line)

        if match:
            operation = match.group(1).strip()
            current["times"][operation] = float(match.group(2))

        # Store completed command.
        if line.startswith("Exit code:"):

            exit_code = int(line.split(":")[1].strip())

            if exit_code != 0:
                successful = False
                print("WARNING: Failed benchmark:", current)

            elif current["target"] == "perf":

                key = (
                    current["N"],
                    current["L"],
                    current["OPT"],
                    current["OMP"]
                )

                perf_data[key] = current["times"]

            else:

                key = (current["N"], current["L"])
                openfhe_data[key] = current["times"]

            current = None


if not successful:
    raise RuntimeError("Some benchmark commands failed. Check the log.")

# --------------------------------------------------
# Operation mappings
# --------------------------------------------------

PRIMITIVES = {
    "CBD + Add": "CBD + Add",
    "Tern. Samp.": "Ternary Sampler",
    "Elem. Wise": "Element-Wise",
    "FFT": "FFT",
    "IFFT": "IFFT",
    "NTT": "NTT",
    "INTT": "INTT"
}

OPERATIONS = {
    "Context Gen.": "ContextGen",
    "Key Gen.": "KeyGen",
    "Relin. Key Gen.": "Relinearization KeyGen",
    "Conjug. Key Gen.": "Conjugation KeyGen",
    "Inner. Key Gen.": "Inner Product KeyGen",
    "Encode": "Encode",
    "Decode": "Decode",
    "Encrypt": "Encrypt",
    "Decrypt": "Decrypt"
}

OPENFHE_OPERATIONS = {
    "Context Gen.": ("ContextGen", "Context Generation"),
    "Key Gen.": ("KeyGen", "KeyGen"),
    "Relin. Key Gen.": (
        "Relinearization KeyGen",
        "Relinearization KeyGen"
    ),
    "Conjug. Key Gen.": (
        "Conjugation KeyGen",
        "Conjugation KeyGen"
    ),
    "Inner. Key Gen.": (
        "Inner Product KeyGen",
        "Inner Product KeyGen"
    ),
    "Encode": ("Encode", "Encode"),
    "Encrypt": ("Encrypt", "Encrypt")
}


def get_time(data, key, operation):
    """Retrieve an operation's execution time."""

    if key not in data:
        raise KeyError(f"Missing benchmark configuration: {key}")

    if operation not in data[key]:
        raise KeyError(
            f"Missing operation '{operation}' for configuration {key}"
        )

    return data[key][operation]


def speedup(baseline, optimized):
    """Calculate speedup."""

    if optimized <= 0:
        raise ValueError("Execution time must be positive.")

    return baseline / optimized


def fmt(value):
    return f"{value:.2f}"


def check_completeness():
    """Verify that all expected configurations exist."""

    for N in N_VALUES:
        for L in L_VALUES:

            for OPT, omp_values in [
                (0, [1]),
                (1, THREADS)
            ]:

                for OMP in omp_values:
                    key = (N, L, OPT, OMP)

                    if key not in perf_data:
                        raise KeyError(
                            f"Missing performance test: {key}"
                        )

            key = (N, L)

            if key not in openfhe_data:
                raise KeyError(
                    f"Missing OpenFHE test: {key}"
                )


# --------------------------------------------------
# LaTeX helper functions
# --------------------------------------------------

def operation_header(name):
    return (
        r"\multirow{3}{*}{\begin{tabular}[c]{@{}l@{}}"
        rf"\textbf{{{name}}}"
        r"\end{tabular}}"
    )


def write_thread_rows(name, values, grouped=False):
    """Generate three rows, optionally inserting spacer cells between L groups."""
    result = []
    for i, thread in enumerate(THREADS):
        prefix = operation_header(name) if i == 0 else ""
        thread_values = values[thread]

        if grouped:
            groups = [
                thread_values[j:j + len(N_VALUES)]
                for j in range(0, len(thread_values), len(N_VALUES))
            ]
            values_str = " &  & ".join(
                " & ".join(f"${fmt(v)}$" for v in group)
                for group in groups
            )
        else:
            values_str = " & ".join(f"${fmt(v)}$" for v in thread_values)

        result.append(f"{prefix}\n & ${thread}$ & {values_str} " + r"\\")

    return "\n".join(result)


def write_group_rows(name, values):
    """Generate one OpenFHE comparison row with spacer cells between L groups."""
    groups = [
        values[j:j + len(N_VALUES)]
        for j in range(0, len(values), len(N_VALUES))
    ]
    values_str = " &  & ".join(
        " & ".join(f"${fmt(v)}$" for v in group)
        for group in groups
    )
    return (
        r"\multirow{1}{*}{\begin{tabular}[c]{@{}l@{}}"
        rf"\textbf{{{name}}}"
        r"\end{tabular}}"
        f" & {values_str} " + r"\\"
    )


# --------------------------------------------------
# TABLE 1
# Primitive speedups
# --------------------------------------------------

def generate_table1():

    rows = []

    for label, operation in PRIMITIVES.items():

        values = {}

        for thread in THREADS:

            values[thread] = []

            for N in N_VALUES:

                baseline_key = (N, PRIMITIVE_L, 0, 1)
                optimized_key = (N, PRIMITIVE_L, 1, thread)

                baseline = get_time(
                    perf_data, baseline_key, operation
                )

                optimized = get_time(
                    perf_data, optimized_key, operation
                )

                values[thread].append(
                    speedup(baseline, optimized)
                )

        rows.append(write_thread_rows(label, values))
        rows.append(r"\midrule")

    rows.pop()

    return r"""
\begin{table*}[!ht]
\caption{Speedup of the optimized primitives over the baseline for a variable number of threads.}
\label{tab::evaluation::speedups_prim}
\centering
\definecolor{vscleanColor}{rgb}{0.3, 0.3, 0.3}
\providecommand{\vsClean}[1]{\color{vscleanColor}\small ($#1$)}
\begin{tabular}{lrrrr}
\toprule
\multirow{2}{*}{\textbf{Operations}} &
\multirow{2}{*}{${\textbf{Threads}}$} &
\multicolumn{3}{c}{\textbf{N}} \\
\cline{3-5}
&& $2^{13}$ & $2^{15}$ & $2^{17}$ \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table*}
"""


# --------------------------------------------------
# TABLE 2
# CKKS operation speedups
# --------------------------------------------------

def generate_table2():

    rows = []

    for label, operation in OPERATIONS.items():

        values = {}

        for thread in THREADS:

            values[thread] = []

            for L in L_VALUES:

                for N in N_VALUES:

                    baseline_key = (N, L, 0, 1)
                    optimized_key = (N, L, 1, thread)

                    baseline = get_time(
                        perf_data, baseline_key, operation
                    )

                    optimized = get_time(
                        perf_data, optimized_key, operation
                    )

                    values[thread].append(
                        speedup(baseline, optimized)
                    )

        # Each thread has nine values: three N values per L.
        rows.append(write_thread_rows(label, values, grouped=True))
        rows.append(r"\midrule")

    rows.pop()

    # Construct the three L headers.
    header = (
        r"\begin{tabular}{lrrrrc@{~~}rrrc@{~~}rrr}"
        "\n"
        r"\toprule"
        "\n"
        r"\multirow{2}{*}{\textbf{Operations}} & "
        r"\multirow{2}{*}{\textbf{Threads}} & "
    )

    header += " & & ".join(
        rf"\multicolumn{{3}}{{c}}{{\textbf{{N (L={L})}}}}"
        for L in L_VALUES
    )

    header += r" \\" + "\n"

    header += (
        r"\cline{3-5}\cline{7-9}\cline{11-13}"
        "\n&& "
    )

    header += (
        r" & ".join(
            [r"$2^{13}$", r"$2^{15}$", r"$2^{17}$"]
            for _ in []
        ) if False else
        " & & ".join(
            [r"$2^{13}$ & $2^{15}$ & $2^{17}$"] * 3
        )
    )

    header += r" \\" + "\n\midrule" + "\n"

    return r"""
\begin{table*}[!ht]
\caption{Speedup of the optimized CKKS operations over the baseline for a variable number of threads.}
\label{tab::evaluation::speedups_ops}
\centering
\definecolor{vscleanColor}{rgb}{0.3, 0.3, 0.3}
\providecommand{\vsClean}[1]{\color{vscleanColor}\small ($#1$)}
\resizebox{\linewidth}{!}{
""" + header + "\n".join(rows) + r"""
\bottomrule
\end{tabular}}
\end{table*}
"""


# --------------------------------------------------
# TABLE 3
# Comparison against OpenFHE
# --------------------------------------------------

def generate_table3():

    rows = []

    for label, (perf_operation, fhe_operation) in OPENFHE_OPERATIONS.items():

        values = []

        for L in L_VALUES:

            for N in N_VALUES:

                # Optimized implementation, four threads.
                perf_key = (
                    N, L, 1, OPENFHE_COMPARE_THREADS
                )

                # OpenFHE configuration.
                fhe_key = (N, L)

                optimized = get_time(
                    perf_data, perf_key, perf_operation
                )

                openfhe = get_time(
                    openfhe_data, fhe_key, fhe_operation
                )

                values.append(speedup(openfhe, optimized))

        rows.append(write_group_rows(label, values))
        rows.append(r"\midrule")

    # Combine Decrypt and Decode for a fair comparison.
    decrypt_values = []

    for L in L_VALUES:
        for N in N_VALUES:

            perf_key = (
                N, L, 1, OPENFHE_COMPARE_THREADS
            )

            fhe_key = (N, L)

            decrypt = get_time(
                perf_data, perf_key, "Decrypt"
            )

            decode = get_time(
                perf_data, perf_key, "Decode"
            )

            openfhe = get_time(
                openfhe_data, fhe_key, "Decrypt + Decode"
            )

            decrypt_values.append(
                speedup(openfhe, decrypt + decode)
            )

    rows.append(
        write_group_rows("Decrypt + Decode", decrypt_values)
    )

    header = (
        r"\begin{tabular}{lrrrc@{~~}rrrc@{~~}rrr}"
        "\n"
        r"\toprule"
        "\n"
        r"\multirow{2}{*}{\textbf{Operations}} & "
    )

    header += " & & ".join(
        rf"\multicolumn{{3}}{{c}}{{\textbf{{N (L={L})}}}}"
        for L in L_VALUES
    )

    header += r" \\" + "\n"

    header += (
        r"\cline{2-4}\cline{6-8}\cline{10-12}"
        "\n& "
        + " & & ".join(
            [r"$2^{13}$ & $2^{15}$ & $2^{17}$"] * 3
        )
        + r" \\" + "\n\midrule"
    )

    return r"""
\begin{table*}[!ht]
\caption{Speedup of the optimized CKKS operations over OpenFHE for $4$ threads.}
\label{tab::evaluation::speedups_openfhe}
\centering
\definecolor{vscleanColor}{rgb}{0.3, 0.3, 0.3}
\providecommand{\vsClean}[1]{\color{vscleanColor}\small ($#1$)}
\resizebox{\linewidth}{!}{
""" + header + "\n" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}}
\end{table*}
"""


# --------------------------------------------------
# Main
# --------------------------------------------------

check_completeness()

table1 = generate_table1()
table2 = generate_table2()
table3 = generate_table3()

with open(OUTPUT_FILE, "w") as output:

    output.write(table1)
    output.write("\n\n")
    output.write(table2)
    output.write("\n\n")
    output.write(table3)

print("LaTeX tables successfully generated!")
print(f"Output file: {OUTPUT_FILE}")
print(f"Performance configurations: {len(perf_data)}")
print(f"OpenFHE configurations: {len(openfhe_data)}")
print(f"Primitive table L: {PRIMITIVE_L}")
print(f"OpenFHE comparison threads: {OPENFHE_COMPARE_THREADS}")