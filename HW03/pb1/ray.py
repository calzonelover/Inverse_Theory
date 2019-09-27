# Reference
# https://www.doc.ic.ac.uk/~dfg/graphics/graphics2008/GraphicsLecture09.pdf

def ray(nx, ny, dx, x1, z1, x2, z2):
    # case 1 same
    if x1 == x2 and z1 == z2:
        pass
    # case 2 ver
    elif x1 == x2:
        pass
    # case 3 hor
    elif z1 == z2:
        pass
    # case 4 tilt
    else:
        m = (z2-z1)/(x2-x1)
        c = z1 - m*x1
        

# def path(s_i, r_i, R):
#     pass # return vector range

# def 