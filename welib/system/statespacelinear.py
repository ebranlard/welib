import numpy as np

# --------------------------------------------------------------------------------}
# --- Simple statespace functions
# --------------------------------------------------------------------------------{
def lti_state_space_function(t, x, u, p):
    return p['A'].dot(x) + p['B'].dot(u)

def lti_output_function(t, x, u, p):
    return p['C'].dot(x) + p['D'].dot(u)
