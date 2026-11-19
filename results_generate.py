import subprocess
from datetime import datetime

# Configuration
N_VALUES = [2**13, 2**15, 2**17]
L_VALUES = [8, 16, 24]

OMP_VALUES = [1, 2, 4, 8]
OPT_VALUES = [0, 1]

N_TESTS = 1

LOG_FILE = "benchmark_results.log"


def run_command(command, log_file):
    """Execute a command and redirect its output to the log file."""

    command_str = " ".join(command)

    log_file.write("\n" + "=" * 80 + "\n")
    log_file.write(f"Command: {command_str}\n")
    log_file.write(f"Started at: {datetime.now()}\n")
    log_file.write("=" * 80 + "\n")
    log_file.flush()

    print(f"\nRunning: {command_str}")

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1
    )

    for line in process.stdout:
        print(line, end="")
        log_file.write(line)
        log_file.flush()

    process.wait()

    log_file.write(f"\nExit code: {process.returncode}\n")
    log_file.write(f"Finished at: {datetime.now()}\n")
    log_file.flush()

    if process.returncode != 0:
        print(f"WARNING: Command failed with exit code {process.returncode}")


with open(LOG_FILE, "w") as log:

    log.write(f"Benchmark execution started at {datetime.now()}\n")

    # --------------------------------------------------
    # 1. Performance tests
    # --------------------------------------------------

    log.write("\n########## PERFORMANCE TESTS ##########\n")

    for N in N_VALUES:
        N2 = N // 2

        for L in L_VALUES:

            for OPT in OPT_VALUES:

                # When OPT=0, OMP_NUM_THREADS is always 1.
                omp_values = [1] if OPT == 0 else OMP_VALUES

                for OMP in omp_values:

                    command = [
                        "make",
                        "TEST=perf",
                        "ARCH=arm",
                        f"OPT={OPT}",
                        f"OMP_NUM_THREADS={OMP}",
                        f"N_TESTS={N_TESTS}",
                        f"N={N}",
                        f"N2={N2}",
                        f"L={L}"
                    ]

                    run_command(command, log)

    # --------------------------------------------------
    # 2. OpenFHE tests
    # --------------------------------------------------

    log.write("\n########## OPENFHE TESTS ##########\n")

    for N in N_VALUES:
        N2 = N // 2

        for L in L_VALUES:

            command = [
                "make",
                "openfhe",
                f"N={N}",
                f"N2={N2}",
                f"L={L}",
                "OMP_NUM_THREADS=4",
                f"N_TESTS={N_TESTS}"
            ]

            run_command(command, log)

    log.write(f"\nAll benchmarks finished at {datetime.now()}\n")

print(f"\nAll benchmarks finished. Results saved to {LOG_FILE}")
