from fun2.logging import log_step
import subprocess

# -----------------------------------------------------------------------------
# Helper Function: Execute Shell Commands
# -----------------------------------------------------------------------------
@log_step
def run_single_command(command: str, working_dir: str = ".") -> bytes:
    """Executes a shell command and returns its stdout as bytes.

    Args:
        command (str): The shell command to execute.
        working_dir (str, optional): The directory to execute the command in. Defaults to '.'.

    Returns:
        bytes: The standard output from the command.

    Raises:
        RuntimeError: If the command exits with a non-zero status.
    """
    process = subprocess.Popen(
        command,
        shell=True,
        cwd=working_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout_data, stderr_data = process.communicate()
    if process.returncode != 0:
        error_msg = stderr_data.decode("utf-8", errors="ignore").strip()
        raise RuntimeError(
            f"Command failed (exit code {process.returncode}): {command}\n{error_msg}"
        )
    return stdout_data