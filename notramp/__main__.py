import os, sys
import notramp.notramp_main

if __name__ == "__main__":
    sys.path.append(os.path.split(__file__)[0])
    notramp.notramp_main.run_notramp()
