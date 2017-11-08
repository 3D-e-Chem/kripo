from math import sqrt, pi, cos, sin


def distance_between_positions(pos1, pos2):
    d = (
        pos1[0] - pos2[0],
        pos1[1] - pos2[1],
        pos1[2] - pos2[2],
    )
    return sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])


def cross_product(a, b):
    """returns the cross product of two vectors"""
    return [a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]]


def center_of_triangle(atom1, atom2, atom3):
    pos1 = atom1.location()
    pos2 = atom2.location()
    pos3 = atom3.location()
    a12 = (
        pos1[0] - pos2[0],
        pos1[1] - pos2[1],
        pos1[2] - pos2[2],
    )
    a23 = (
        pos2[0] - pos3[0],
        pos2[1] - pos3[1],
        pos2[2] - pos3[2],
    )
    d = cross_product(a12, a23)
    d_dist = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
    return (
        d[0] / d_dist,
        d[1] / d_dist,
        d[2] / d_dist,
    )


def vector_rotate(v, angle, axis):
    """Rotate the vector v by angle (in radians) about the axis.
    """
    # NORMALIZE
    d = 0
    for i in range(3):
        d = d + (v[i] * v[i])
    d = sqrt(d)
    v = [x/d for x in v]

    # CONVERT ANGLE TO RADIANS
    angle = (float(angle)/360.0)*(2*pi)

    # GET SINE AND COSINE
    c = cos(angle)
    s = sin(angle)

    # AXIS
    (n0, n1, n2) = axis

    # ROTATION MATRIX
    r = [
        [c + (1-c)*n0*n0, (1-c)*n0*n1 - s*n2, (1-c)*n0*n2 + s*n1],
        [(1-c)*n1*n0 + s*n2, c + (1-c)*n1*n1, (1-c)*n1*n2 - s*n0],
        [(1-c)*n2*n0 - s*n1, (1-c)*n2*n1 + s*n0, c + (1-c)*n2*n2]
    ]

    # TRANSFORM VECTOR
    v2 = []
    for i in range(3):
        v2.append(0)
        for j in range(3):
            v2[i] += r[i][j] * v[j]

    return v2
