import numpy as np
# import matplotlib.pyplot as plt

# import pdb


def trapz_cond(x, y, c=None):
    """
    Do trapz over the function y,
    possibly conditioned upon condition c (array of bool).

    All must be equal-sized arrays
    """
    # Make sure things are arrays
    x, y = np.array(x), np.array(y)
    if c is not None:
        assert type(c) is np.ndarray # must be array

        if c.dtype == np.bool_:
            c = c.astype(float) - 0.5
        elif np.issubdtype(c.dtype, np.number):
            c = np.array(c, dtype=float)
        else:
            raise TypeError('Invalid datatype for condition array c in trapz_cond.')
    else:
        c = np.ones(len(x))
    
    # Arrange things in proper order
    dx = np.diff(x)
    if (dx < 0).all():
        x, y, c = np.flip(x), np.flip(y), np.flip(c)
        dx = -np.flip(dx)
    elif not (dx > 0).all():
        raise ValueError('x must be monotonically increasing or decreasing.')

    # Do conditional integral (L, R := left, right)
    y_dx = np.zeros(len(dx))
    yL, yR = y[:-1], y[1:]
    cL, cR = c[:-1], c[1:]
    sgn_cL, sgn_cR = np.sign(cL), np.sign(cR)

    cond = (sgn_cL >= 0) & (sgn_cR >= 0) # both L & R satisfy condition (or it's marginal)
    y_dx[cond] = 0.5 * (yL[cond] + yR[cond]) * dx[cond]
    
    cond = (sgn_cL > 0) & (sgn_cR < 0) # L satisfies condition, but not R
    dxc = cL[cond] * dx[cond] / (cL[cond] - cR[cond])
    yc = yL[cond] * (1 - dxc / dx[cond]) + yR[cond] * dxc / dx[cond]
    y_dx[cond] = 0.5 * (yL[cond] + yc) * dxc

    cond = (sgn_cL < 0) & (sgn_cR > 0) # R satisfies condition, but not L
    dxc = cL[cond] * dx[cond] / (cL[cond] - cR[cond])
    yc = yL[cond] * (1 - dxc / dx[cond]) + yR[cond] * dxc / dx[cond]
    y_dx[cond] = 0.5 * (yc + yR[cond]) * (dx[cond] - dxc)

    int_y_dx = np.sum(y_dx)

    return int_y_dx



# x = np.linspace(0, 5, 10)
# y = x ** 2
# c = x - 2

# x = np.array([0, 1])
# y = np.array([0, 1])
# c = np.array([1, -1])

# int_y_dx_old = []
# int_y_dx_new = []
# int_y_dx_newf = []
# for ii in range(len(x)):
#     int_y_dx_old.append(np.trapz(x=x[:ii+1][c[:ii+1]>0], y=y[:ii+1][c[:ii+1]>0]))
#     int_y_dx_new.append(cumtrapz_cond(x[:ii+1], y[:ii+1], c[:ii+1]))

#     int_y_dx_newf.append(cumtrapz_cond(np.flip(x[:ii+1]), np.flip(y[:ii+1]), np.flip(c[:ii+1])))
# # int_y_dx_the = x ** 3 / 3 - 8 / 3
# # int_y_dx_the[c<0] = 0

# plt.close()
# plt.scatter(x, int_y_dx_old, label='old', s=1)
# plt.scatter(x, int_y_dx_new, label='new', s=1)
# plt.scatter(x, int_y_dx_newf, label='newf', s=1)
# # plt.scatter(x, int_y_dx_the, label='theory', s=1)
# plt.legend(loc='upper left')
# plt.show()






















