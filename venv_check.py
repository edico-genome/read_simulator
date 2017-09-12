import sys

def is_venv():
    """ confirm we are running in a virtualenv """
    return (hasattr(sys, 'real_prefix') or
            (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix))

if not is_venv():
    print "Please run in appropriate Python2.7 venv"
    sys.exit(1)
