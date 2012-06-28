from amcmorl_py_tools.proc_tools import arg_ind_fn_where

def test_arg_ind_fn_where():
    a = np.array([16, 43, 96, 68, 84, 93, 46, 78, 12, 16]);
    cond = a < 50
    np.testing.assert_equal(arg_ind_fn_where(np.argmax, a, cond), 6)
