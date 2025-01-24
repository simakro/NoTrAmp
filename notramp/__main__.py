import os, sys
import notramp.notramp_main

if __name__ == "__main__":
    pkg_path = os.path.split(__file__)[0]
    if pkg_path not in sys.path:
        sys.path.append(pkg_path)
    notramp.notramp_main.run_notramp()
