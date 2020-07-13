import atexit

def pwscf_finalise():
    pwpy_pwscf_finalise()


atexit.register(pwscf_finalise)
