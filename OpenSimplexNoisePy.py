import random

__name__ = "OpenSimplexNoisePy"

# Stretch constant to skew coords onto lattice
STRETCH_CONSTANT_3D = -1.0 / 6

#Squish constant to skew coords off of lattice
SQUISH_CONSTANT_3D = 1.0 / 3

#testcomment

#Default Perm Table
PERM_DEFAULT = [151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225,
                140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148,
                247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32,
                57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175,
                74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122,
                60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54,
                65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169,
                200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64,
                52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212,
                207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213,
                119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
                129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104,
                218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241,
                81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157,
                184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93,
                222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180]

#Table of 3D gradient values
#Updated on 9/19/14
gradients3D = [0, 3, 2, 0, 2, 3, 3, 0, 2, 2, 0, 3, 3, 2, 0, 2, 3, 0,
               0, -3, 2, 0, 2, -3, -3, 0, 2, 2, 0, -3, -3, 2, 0, 2, -3, 0,
               0, 3, -2, 0, -2, 3, 3, 0, -2, -2, 0, 3, 3, -2, 0, -2, 3, 0,
               0, -3, -2, 0, -2, -3, -3, 0, -2, -2, 0, -3, -3, -2, 0, -2, -3, 0]

#Finish declaring global variables
permGradIndex3D = [0 for x in range(256)]
perm = [0 for x in range(256)]


#Set up the perm tables and the gradient index table
def setUp(default, seed):
    global perm
    global permGradIndex3D
    global PERM_DEFAULT
    testBool = True

    if (default):
        perm = PERM_DEFAULT
        for i in range(256):
            permGradIndex3D[i] = (int)((perm[i] % (len(gradients3D) / 3)) * 3)
    else:
        source = [0 for x in range(256)]
        random.seed(seed)

        for i in range(256):
            source[i] = i

        for i in range(255, -1, -1):
            r = random.randint(0, i + 1)
            perm[i] = source[r]
            permGradIndex3D[i] = (int)((perm[i] % (len(gradients3D) / 3)) * 3)
            source[r] = source[i]


def fastFloor(x):
    xi = int(x)
    return (xi - 1) if (x < xi) else (xi)


def extrapolate(xsb, ysb, zsb, dx, dy, dz):
    global permGradIndex3D
    global gradients3D
    global perm

    index = permGradIndex3D[(perm[(perm[xsb & 0xFF] + ysb) & 0xFF] + zsb) & 0xFF]
    print((gradients3D[index] * dx) + (gradients3D[index + 1] * dy) + (gradients3D[index + 2] * dz))
    return (gradients3D[index] * dx) + (gradients3D[index + 1] * dy) + (gradients3D[index + 2] * dz)


def eval(x, y, z):
    global STRETCH_CONSTANT_3D
    global SQUISH_CONSTANT_3D

    #Skew coords onto simplectic lattice
    stretchOffset = (x + y + z) * STRETCH_CONSTANT_3D
    xs = x + stretchOffset
    ys = y + stretchOffset
    zs = z + stretchOffset

    #Floor to get lattice coords of rhombohedron
    xsb = fastFloor(xs)
    ysb = fastFloor(ys)
    zsb = fastFloor(zs)

    #Squish coords back to rhombo coords
    squishOffset = (x + y + z) * SQUISH_CONSTANT_3D
    xb = x + squishOffset
    yb = y + squishOffset
    zb = z + squishOffset

    #Get lattice coords relative to rhombohedral origin
    xins = xs - xsb
    yins = ys - ysb
    zins = zs - zsb

    #Sum relative coords to figure out what rhombohedron we're inside
    inSum = xins + yins + zins

    #Figure out relative coords to origin
    dx0 = x - xb
    dy0 = y - yb
    dz0 = z - zb

    #reset return value
    value = 0

    if (inSum <= 1):
        aPoint = 0x01
        aScore = xins

        bPoint = 0x02
        bScore = yins

        if (aScore >= bScore and zins > bScore):
            bScore = zins
            bPoint = 0x04
        elif (aScore < bScore and zins > aScore):
            aScore = zins
            aPoint = 0x04

        #Determine the two lattice points not part of the tetrahedron that are closest
        #Depends on closest two tetrahedral vertices
        wins = 1 - inSum
        if (wins > aScore or wins > bScore):
            c = (bPoint) if (bScore > aScore) else (aPoint)

            if ((c & 0x01) == 0):
                xsv_ext0 = xsb - 1
                xsv_ext0 = xsb
                dx_ext0 = dx0 + 1
                dx_ext1 = dx0
            else:
                xsv_ext0 = xsv_ext1 = xsb + 1
                dx_ext0 = dx_ext1 = dx0 - 1

            if ((c & 0x02) == 0):
                ysv_ext0 = ysv_ext1 = ysb
                dy_ext0 = dy_ext1 = dy0

                if ((c & 0x01) == 0):
                    ysv_ext1 -= 1
                    dy_ext1 += 1
                else:
                    ysv_ext0 -= 1
                    dy_ext0 += 1
            else:
                ysv_ext0 = ysv_ext1 = ysb + 1
                dy_ext0 = dy_ext1 = dy0 - 1

            if ((c & 0x04) == 0):
                zsv_ext0 = zsb
                zsv_ext1 = zsb - 1
                dz_ext0 = dz0
                dz_ext1 = dz0 + 1
            else:
                zsv_ext0 = zsv_ext1 = zsb + 1
                dz_ext0 = dz_ext1 = dz0 - 1

        else:  # (0,0,0) is not one of the closest two vertices
            c = (aPoint | bPoint)

            if ((c & 0x01) == 0):
                xsv_ext0 = xsb
                xsv_ext1 = xsb - 1
                dx_ext0 = dx0 - 2 * SQUISH_CONSTANT_3D
                dx_ext1 = dx0 + 1 - SQUISH_CONSTANT_3D
            else:
                xsv_ext0 = xsv_ext1 = xsb + 1
                dx_ext0 = dx0 - 1 - 2 * SQUISH_CONSTANT_3D
                dx_ext1 = dx0 - 1 - SQUISH_CONSTANT_3D

            if ((c & 0x02) == 0):
                ysv_ext0 = ysb
                ysv_ext1 = ysb - 1
                dy_ext0 = dy0 - 2 * SQUISH_CONSTANT_3D
                dy_ext1 = dy0 + 1 - SQUISH_CONSTANT_3D
            else:
                ysv_ext0 = ysv_ext1 = ysb + 1
                dy_ext0 = dy0 - 1 - 2 * SQUISH_CONSTANT_3D
                dy_ext1 = dy0 - 1 - SQUISH_CONSTANT_3D

            if ((c & 0x04) == 0):
                zsv_ext0 = zsb
                zsv_ext1 = zsb - 1
                dz_ext0 = dz0 - 2 * SQUISH_CONSTANT_3D
                dz_ext1 = dz0 + 1 - SQUISH_CONSTANT_3D
            else:
                zsv_ext0 = zsv_ext1 = zsb + 1
                dz_ext0 = dz0 - 1 - 2 * SQUISH_CONSTANT_3D
                dz_ext1 = dz0 - 1 - SQUISH_CONSTANT_3D

        #Contribution (0,0,0)
        attn0 = 2 - dx0 * dx0 - dy0 * dy0 - dz0 * dz0
        if (attn0 > 0):
            attn0 *= attn0
            value = attn0 * attn0 * extrapolate(xsb + 0, ysb + 0, zsb + 0, dx0, dy0, dz0)

        #Contribution (0,0,1)
        dx1 = dx0 - 1 - SQUISH_CONSTANT_3D
        dy1 = dy0 - 0 - SQUISH_CONSTANT_3D
        dz1 = dz0 - 0 - SQUISH_CONSTANT_3D
        attn1 = 2 - dx1 * dx1 - dy1 * dy1 - dz1 * dz1
        if (attn1 > 0):
            attn1 *= attn1
            value += attn1 * attn1 * extrapolate(xsb + 1, ysb + 0, zsb + 0, dx1, dy1, dz1)

        #Contribution (0,1,0)
        dx2 = dx0 - 0 - SQUISH_CONSTANT_3D
        dy2 = dy0 - 1 - SQUISH_CONSTANT_3D
        dz2 = dz1
        attn2 = 2 - dx2 * dx2 - dy2 * dy2 - dz2 * dz2
        if (attn2 > 0):
            attn2 *= attn2
            value += attn2 * attn2 * extrapolate(xsb + 0, ysb + 1, zsb + 0, dx2, dy2, dz2)

        #Contribution (1,0,0)
        dx3 = dx2
        dy3 = dy1
        dz3 = dz0 - 1 - SQUISH_CONSTANT_3D
        attn3 = 2 - dx3 * dx3 - dy3 * dy3 - dz3 * dz3
        if (attn3 > 0):
            attn3 *= attn3
            value += attn3 * attn3 * extrapolate(xsb + 0, ysb + 0, zsb + 1, dx3, dy3, dz3)

    elif (inSum >= 2):
        value = 400

    return value / 28.25


setUp(True, 0)

for y in range(256):
    for x in range(256):
        (eval(x / 256, y / 256, 0))








