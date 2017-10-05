import sys

""" confirm that we are running in a virtualenv """
test1 = hasattr(sys, 'real_prefix')
test2 = hasattr(sys, 'base_prefix') and (sys.base_prefix != sys.prefix)
if not (test1 or test2):
    print "Please run in appropriate Python2.7 venv"
    sys.exit(1)
