import os
import sys
import subprocess


def main():
    # Resolve package root and scooti.sh path
    here = os.path.dirname(os.path.abspath(__file__))
    script = os.path.join(here, 'scooti.sh')
    if not os.path.exists(script):
        sys.stderr.write(f"scooti.sh not found at {script}\n")
        return 1
    # Ensure executable
    try:
        os.chmod(script, 0o755)
    except Exception:
        pass
    # Delegate to bash script with same args
    cmd = ['bash', script] + sys.argv[1:]
    return subprocess.call(cmd)


if __name__ == '__main__':
    raise SystemExit(main())

