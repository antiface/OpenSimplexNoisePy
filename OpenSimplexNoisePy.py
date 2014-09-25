import random
import math

STRETCH_CONSTANT_3D = -1.0 / 6
SQUISH_CONSTANT_3D = 1.0 /3

gradients3D = [0,1,2,
				0,2,1,
				1,0,2,
				2,0,1,
				2,1,0,
				0,-1,2,
				0,2,-1,
				1,0,2,
				2,0,-1,
				-1,2,0,
				2,-1,0,
				0,1,-2,
				0,-2,1,
				1,0,-2,
				-2,0,1,
				1,-2,0,
				-2,1,0,
				0,-1,-2,
				0,-2,-1,
				-1,0,-2,
				-2,0,-1,
				-1,-2,0,
				-2,-1,0]

perm_default = [151,160,137, 91, 90, 15,131, 13,201, 95, 96, 53,194,233,  7,225,
	140, 36,103, 30, 69,142,  8, 99, 37,240, 21, 10, 23,190,  6,148,
	247,120,234, 75,  0, 26,197, 62, 94,252,219,203,117, 35, 11, 32,
	 57,177, 33, 88,237,149, 56, 87,174, 20,125,136,171,168, 68,175,
	 74,165, 71,134,139, 48, 27,166, 77,146,158,231, 83,111,229,122,
	 60,211,133,230,220,105, 92, 41, 55, 46,245, 40,244,102,143, 54,
	 65, 25, 63,161,  1,216, 80, 73,209, 76,132,187,208, 89, 18,169,
	200,196,135,130,116,188,159, 86,164,100,109,198,173,186,  3, 64,
	 52,217,226,250,124,123,  5,202, 38,147,118,126,255, 82, 85,212,
	207,206, 59,227, 47, 16, 58, 17,182,189, 28, 42,223,183,170,213,
	119,248,152,  2, 44,154,163, 70,221,153,101,155,167, 43,172,  9,
	129, 22, 39,253, 19, 98,108,110, 79,113,224,232,178,185,112,104,
	218,246, 97,228,251, 34,242,193,238,210,144, 12,191,179,162,241,
	 81, 51,145,235,249, 14,239,107, 49,192,214, 31,181,199,106,157,
	184, 84,204,176,115,121, 50, 45,127,  4,150,254,138,236,205, 93,
	222,114, 67, 29, 24, 72,243,141,128,195, 78, 66,215, 61,156,180]

perm = [0 for x in range(256)]
permGradIndex3D = [0 for x in range(256)]
source = [0 for x in range(256)]

for i in range(0, 256):
		source[i] = i

random.seed()

for i in range(255, -1, -1):
	r = random.randrange(i+1)
	perm[i] = source[r]
	permGradIndex3D[i] = ((perm[i] % (len(gradients3D) / 3)) * 3)
	source[r] = source[i]

def extrapolate(xsb, ysb, zsb, dx, dy, dz):
	print("xsb:", xsb, ", ysb:", ysb, ", zsb:", zsb);
	index = (perm[(perm[(perm[xsb & 0xFF] + ysb) & 0xFF] + zsb) & 0xFF] % 69)
	#if (index + 2 < (len(gradients3D) - 1)):
	print(index)
	print(len(gradients3D))
	return gradients3D[index] * dx + gradients3D[index + 1] * dy + gradients3D[index + 2] * dz
	#else:
	#	index = (len(gradients3D) - 3)
	#	return gradients3D[index] * dx + gradients3D[index + 1] * dy + gradients3D[index + 2] * dz


def fastFloor(x):
	xi = int(x)
	return (xi - 1) if (x < xi) else xi

def eval(x, y, z):
	STRETCH_CONSTANT_3D = (-1.0 / 6)
	SQUISH_CONSTANT_3D = (1.0 /3)

	#Place coords on a simplectic lattice
	strechOffset = (x + y + z) * STRETCH_CONSTANT_3D
	xs = x + strechOffset
	ys = y + strechOffset
	zs = z + strechOffset

	#Floor to get coords of rhombohedron
	xsb = fastFloor(xs)
	ysb = fastFloor(ys)
	zsb = fastFloor(zs)

	#Skew out to get actual rhombohedron coords
	squishOffset = (x + y + z) * SQUISH_CONSTANT_3D
	xb = xsb + squishOffset
	yb = ysb + squishOffset
	zb = zsb + squishOffset

	#Compute simplectic lattice coordinates relative to the rhombohedral origin
	xins = xs - ysb
	yins = ys - ysb
	zins = zs - zsb

	#Sum those together to determine cell location
	inSum = xins + yins + zins

	#Positions relative to origin point
	dx0 = x - xb
	dy0 = y - yb
	dz0 = z - zb

	#Defined now for use later
	dx_ext0 = 0
	dy_ext0 = 0
	dz_ext0 = 0

	dx_ext1 = 0
	dy_ext1 = 0
	dz_ext1 = 0

	xsv_ext0 = 0
	ysv_ext0 = 0
	zsv_ext0 = 0

	xsv_ext1 = 0
	ysv_ext1 = 0
	zsv_ext1 = 0
	#finish definition

	value = 0
	if (inSum <= 1): #We're inside of Tetrahedron at (0,0,0)
		
		#Determine which two of (0,0,1), (0,1,0), and (1,0,0) are closest
		aPoint = 0x01
		aScore = xins
		bPoint = 0x02
		bScore = yins
		if(aScore >= bScore and zins > bScore):
			bScore = zins
			bPoint = 0x04
		elif(aScore < bScore and zins > aScore):
			aScore = zins
			aPoint = 0x04

		# Now determine two lattice points not on tetrahedron that may contribute
		wins = 1 - inSum
		if(wins > aScore or wins > bScore): #(0,0,0) is one of the closest tetrahedral vertices
			c = bPoint if (bScore > aScore) else aPoint

			if((c & 0x01) == 0):
				xsv_ext0 = xsb - 1
				xsv_ext1 = xsb
				dx_ext0 = dx0 + 1
				dx_ext1 = dx0
			else:
				xsv_ext0 = xsv_ext1 = xsb + 1
				dx_ext0 = dx_ext1 = dx0 - 1

			if((c & 0x02) == 0):
				ysv_ext0 = ysb - 1
				ysv_ext1 = ysb
				dy_ext0 = dy0 + 1
				dy_ext1 = dy0
			else:
				ysv_ext0 = ysv_ext1 = ysb + 1
				dy_ext0 = dy_ext1 = dy0 - 1

			if((c & 0x04) == 0):
				zsv_ext0 = zsb - 1
				zsv_ext1 = zsb
				dz_ext0 = dz0 + 1
				dz_ext1 = dz0
			else:
				zsv_ext0 = zsv_ext1 = zsb + 1
				dz_ext0 = dz_ext1 = dz0 - 1

		else:
			c = (aPoint | bPoint)

			if((c & 0x01) == 0):
				xsv_ext0 = xsb
				xsv_ext1 = xsb - 1
				dx_ext0 = dx0 - 2 * SQUISH_CONSTANT_3D
				dx_ext1 = dx0 + 1 - SQUISH_CONSTANT_3D
			else:
				xsv_ext0 = xsv_ext1 = xsb + 1
				dx_ext0 = dx0 - 1 - 2 * SQUISH_CONSTANT_3D
				dx_ext1 = dx0 + 1 - SQUISH_CONSTANT_3D

			if((c & 0x02) == 0):
				ysv_ext0 = ysb
				ysv_ext1 = ysb - 1
				dy_ext0 = dy0 - 2 * SQUISH_CONSTANT_3D
				dy_ext1 = dy0 + 1 - SQUISH_CONSTANT_3D
			else:
				ysv_ext0 = ysv_ext1 = ysb + 1
				dy_ext0 = dy0 - 1 - 2 * SQUISH_CONSTANT_3D
				dy_ext1 = dy0 + 1 - SQUISH_CONSTANT_3D

			if((c & 0x04) == 0):
				zsv_ext0 = zsb
				zsv_ext1 = zsb - 1
				dz_ext0 = dz0 - 2 * SQUISH_CONSTANT_3D
				dz_ext1 = dz0 + 1 - SQUISH_CONSTANT_3D
			else:
				zsv_ext0 = zsv_ext1 = zsb + 1
				dz_ext0 = dz0 - 1 - 2 * SQUISH_CONSTANT_3D
				dz_ext1 = dz0 + 1 - SQUISH_CONSTANT_3D

		#contribution (0,0,0)
		attn0 = 2 - (dx0 ** 2) - (dy0 ** 2) - (dz0 ** 2)
		if (attn0 > 0):
			attn0 = attn0 * attn0
			value = attn0 * attn0 * extrapolate(xsb + 0, ysb + 0, zsb + 0, dx0, dy0, dz0)

		#contribution(0,0,1)
		dx1 = dx0 - 1 - SQUISH_CONSTANT_3D
		dy1 = dy0 - 0 - SQUISH_CONSTANT_3D
		dz1 = dz0 - 0 - SQUISH_CONSTANT_3D
		attn1 = 2 - (dx1 ** 2) - (dy1 ** 2) - (dz1 ** 2)
		if (attn1 > 0):
			attn1 = attn1 * attn1
			value = value + (attn1 * attn1) * extrapolate(xsb + 1, ysb + 0, zsb + 0, dx1, dy1, dz1)

		#contribution(0,1,0)
		dx2 = dx0 - 0 - SQUISH_CONSTANT_3D
		dy2 = dy0 - 1 - SQUISH_CONSTANT_3D
		dz2 = dz1
		attn2 = 2 - (dx2 ** 2) - (dy2 ** 2) - (dz2 ** 2)
		if (attn2 > 0):
			attn2 = attn2 * attn2
			value = value + (attn2 * attn2) * extrapolate(xsb + 0, ysb + 1, zsb + 0, dx2, dy2, dz2)

		#contribution(1,0,0)
		dx3 = dx2
		dy3 = dy1
		dz3 = dz0 - 1 - SQUISH_CONSTANT_3D
		attn3 = 2 - (dx3 ** 2) - (dy3 ** 2) - (dz3 ** 2)
		if (attn3 > 0):
			attn3 = attn3 * attn3
			value = value + (attn3 * attn3) * extrapolate(xsb + 0, ysb + 0, zsb + 1, dx3, dy3, dz3)
	elif (inSum >= 2): #We're inside the Tetrahedron at (1,1,1)

		#Determine which vertices are closest out of (1,1,0), (1,0,1), and (0,1,1) but not (1,1,1)
		aPoint = 0x06
		aScore = xins
		bPoint = 0x05
		bScore = yins
		if(aScore <= bScore and zins < bScore):
			bScore = zins
			bPoint = 0x03
		elif(aScore > bScore and zins < aScore):
			aScore = zins
			aPoint = 0x03

		wins = 3 - inSum
		if (wins < aScore or wins < bScore):
			c = bPoint if (bScore < aScore) else aPoint

			if ((c & 0x01) != 0):
				xsv_ext0 = xsb + 2
				xsv_ext1 = xsb + 1
				dx_ext0 = dx0 - 2 - 3 * SQUISH_CONSTANT_3D
				dx_ext1 = dx0 - 1 - 3 * SQUISH_CONSTANT_3D
			else:
				xsv_ext0 = xsv_ext1 = xsb
				dx_ext0 = dx_ext1 = dx0 - 3 * SQUISH_CONSTANT_3D

			if ((c & 0x02) != 0):
				ysv_ext0 = ysv_ext1 = ysb + 1
				dy_ext0 = dy_ext1 = dy0 - 1 - 3 * SQUISH_CONSTANT_3D
				if ((c & 0x01) != 0):
					ysv_ext1 += 1
					dy_ext1 -= 1
				else:
					ysv_ext0 += 1
					dy_ext0 -= 1
			else:
				ysv_ext0 = ysv_ext1 = ysb
				dy_ext0 = dy_ext1 = dy0 - 3 * SQUISH_CONSTANT_3D

			if ((c & 0x04) != 0):
				zsv_ext0 = zsb + 1
				zsv_ext1 = zsb + 2
				dz_ext0 = dz0 - 1 - 3 * SQUISH_CONSTANT_3D
				dz_ext1 = dz0 - 2 - 3 * SQUISH_CONSTANT_3D
			else:
				zsv_ext0 = zsv_ext1 = zsb
				dz_ext0 = dz_ext1 = dz0 - SQUISH_CONSTANT_3D

		else: #(1,1,1 is not one of the closest points)
			c = (aPoint & bPoint)

			if ((c & 0x01) != 0):
				xsv_ext0 = xsb + 1
				xsv_ext1 = xsb + 2
				dx_ext0 = dx0 - 1 - SQUISH_CONSTANT_3D
				dx_ext1 = dx0 - 2 - 2 * SQUISH_CONSTANT_3D
			else:
				xsv_ext0 = xsv_ext1 = xsb
				dx_ext0 = dx0 - SQUISH_CONSTANT_3D
				dx_ext1 = dx0 - 2 * SQUISH_CONSTANT_3D

			if ((c & 0x02) != 0):
				ysv_ext0 = ysb + 1
				ysv_ext1 = ysb + 2
				dy_ext0 = dy0 - 1 - SQUISH_CONSTANT_3D
				dy_ext1 = dy0 - 2 - 2 * SQUISH_CONSTANT_3D
			else:
				ysv_ext0 = ysv_ext1 = ysb
				dy_ext0 = dy0 - SQUISH_CONSTANT_3D
				dy_ext1 = dy0 - 2 * SQUISH_CONSTANT_3D

			if ((c & 0x04) != 0):
				zsv_ext0 = zsb + 1
				zsv_ext1 = zsb + 2
				dz_ext0 = dz0 - 1 - SQUISH_CONSTANT_3D
				dz_ext1 = dz0 - 2 - 2 * SQUISH_CONSTANT_3D
			else:
				zsv_ext0 = zsv_ext1 = ysb
				dz_ext0 = dz0 - SQUISH_CONSTANT_3D
				dz_ext1 = dz0 - 2 * SQUISH_CONSTANT_3D
		
		#Contribution (1,1,0)
		dx3 = dx0 - 1 - 2 * SQUISH_CONSTANT_3D
		dy3 = dy0 - 1 - 2 * SQUISH_CONSTANT_3D
		dz3 = dz0 - 0 - 2 * SQUISH_CONSTANT_3D
		attn3 = 2 - (dx3 ** 2) - (dy3 ** 2) - (dz3 ** 2)
		if (attn3 > 0):
			attn3 = attn3 * attn3
			value = (attn3 * attn3) * extrapolate(xsb + 1, ysb + 1, zsb + 0, dx3, dy3, dz3)

		#Contribution (1,0,1)
		dx2 = dx3
		dy2 = dy0 - 0 - 2 * SQUISH_CONSTANT_3D
		dz2 = dz0 - 1 - 2 * SQUISH_CONSTANT_3D
		attn2 = 2 - (dx2 ** 2) - (dy2 ** 2) - (dz2 ** 2)
		if (attn2 > 0):
			attn2 = attn2 * attn2
			value = value + (attn2 * attn2) * extrapolate(xsb + 1, ysb + 0, zsb + 1, dx2, dy2, dz2)

		#Contribution (0,1,1)
		dx1 = dx0 - 0 - 2 * SQUISH_CONSTANT_3D
		dy1 = dy3
		dz1 = dz2
		attn1 = 2 - (dx1 ** 2) - (dy1 ** 2) - (dz1 ** 2)
		if (attn1 > 0):
			attn1 = attn1 * attn1
			value = value + (attn1 * attn1) * extrapolate(xsb + 0, ysb + 1, zsb + 1, dx1, dy1, dz1) 

		#Contribution (1,1,1)
		dx0 = dx0 - 1 - 3 * SQUISH_CONSTANT_3D
		dy0 = dy0 - 1 - 3 * SQUISH_CONSTANT_3D
		dz0 = dz0 - 1 - 3 * SQUISH_CONSTANT_3D
		attn0 = 2 - (dx0 ** 2) - (dy0 ** 2) - (dz0 ** 2)
		if (attn0 > 0):
			attn0 = attn0 * attn0
			value = value + (attn0 * attn0) * extrapolate(xsb + 1, ysb + 1, zsb + 1, dx0, dy0, dz0)

	else: #We're inside the octahedron in between
		aScore = None
		aPoint = None
		aIsFurtherSide = None
		bScore = None
		bPoint = None
		bIsFurtherSide = None

		p1 = xins + yins
		if (p1 > 1):
			aScore = p1 - 1
			aPoint = 0x03
			aIsFurtherSide = True
		else:
			bScore = 1 - p1
			bPoint = 0x04
			bIsFurtherSide = False

		p2 = xins + zins
		if (p2 > 1):
			bScore = p2 - 1
			bPoint = 0x05
			bIsFurtherSide = True
		else:
			bScore = 1 - p2
			bPoint = 0x02
			bIsFurtherSide = False

		p3 = yins + zins
		if (p3 > 1):
			score = p3 - 1
			if (aScore <= bScore and aScore < score):
				aScore = score
				aPoint = 0x06
				aIsFurtherSide = True
			elif (aScore > bScore and bScore < score):
				bScore = score
				bPoint = 0x06
				bIsFurtherSide = True
		else:
			score = 1 - p3
			if (aScore <= bScore and aScore < score):
				aScore = score
				aPoint = 0x01
				aIsFurtherSide = False
			elif (aScore > bScore and bScore < score):
				bScore = score
				bPoint = 0x01
				bIsFurtherSide = True

		if (aIsFurtherSide == bIsFurtherSide):
			if (aIsFurtherSide == True):
				dx_ext0 = dx0 - 1 - 3 * SQUISH_CONSTANT_3D
				dy_ext0 = dy0 - 1 - 3 * SQUISH_CONSTANT_3D
				dz_ext0 = dz0 - 1 - 3 * SQUISH_CONSTANT_3D
				xsv_ext0 = xsb + 1
				ysv_ext0 = ysb + 1
				zsv_ext0 = zsb + 1

				c = (aPoint & bPoint)

				if ((c & 0x01) != 0):
					dx_ext1 = dx0 - 2 - 2 * SQUISH_CONSTANT_3D
					dy_ext1 = dy0 - 2 * SQUISH_CONSTANT_3D
					dz_ext1 = dz0 - 2 * SQUISH_CONSTANT_3D
					xsv_ext1 = xsb + 2
					ysv_ext1 = ysb
					zsv_ext2 = zsb
				elif ((c & 0x02) != 0):
					dx_ext1 = dx0 - 2 * SQUISH_CONSTANT_3D
					dy_ext1 = dy0 - 2 - 2 * SQUISH_CONSTANT_3D
					dz_ext1 = dz0 - 2 * SQUISH_CONSTANT_3D
					xsv_ext1 = xsb
					ysv_ext1 = ysb + 2
					zsv_ext2 = zsb
				else:
					dx_ext1 = dx0 - 2 * SQUISH_CONSTANT_3D
					dy_ext1 = dy0 - 2 * SQUISH_CONSTANT_3D
					dz_ext1 = dz0 - 2 - 2 * SQUISH_CONSTANT_3D
					xsv_ext1 = xsb
					ysv_ext1 = ysb
					zsv_ext2 = zsb + 2
			else:
				dx_ext0 = dx0
				dy_ext0 = dy0
				dz_ext0 = dz0
				xsv_ext0 = xsb
				ysv_ext0 = ysb
				zsb_ext0 = zsb
				c = (aPoint | bPoint)

				if ((c & 0x01) == 0):
					dx_ext1 = dx0 + 1 - SQUISH_CONSTANT_3D
					dy_ext1 = dy0 - 1 - SQUISH_CONSTANT_3D
					dz_ext1 = dz0 - 1 - SQUISH_CONSTANT_3D
					xsv_ext1 = xsb - 1
					ysv_ext1 = ysb + 1
					zsb_ext1 = zsb + 1
				elif ((c & 0x02) == 0):
					dx_ext1 = dx0 - 1 - SQUISH_CONSTANT_3D
					dy_ext1 = dy0 + 1 - SQUISH_CONSTANT_3D
					dz_ext1 = dz0 - 1 - SQUISH_CONSTANT_3D
					xsv_ext1 = xsb + 1
					ysv_ext1 = ysb - 1
					zsb_ext1 = zsb + 1
				else: 
					dx_ext1 = dx0 - 1 - SQUISH_CONSTANT_3D
					dy_ext1 = dy0 - 1 - SQUISH_CONSTANT_3D
					dz_ext1 = dz0 + 1 - SQUISH_CONSTANT_3D
					xsv_ext1 = xsb + 1
					ysv_ext1 = ysb + 1
					zsb_ext1 = zsb - 1
			
		else:
			c1 = None
			c2 = None
			if(aIsFurtherSide):
				c1 = aPoint
				c2 = bPoint
			else:
				c1 = bPoint
				c2 = aPoint

			if ((c1 & 0x01) == 0):
				dx_ext0 = dx0 + 1 - SQUISH_CONSTANT_3D
				dy_ext0 = dy0 - 1 - SQUISH_CONSTANT_3D
				dz_ext0 = dz0 - 1 - SQUISH_CONSTANT_3D
				xsv_ext0 = xsb - 1
				ysv_ext0 = ysb + 1
				zsb_ext0 = zsb + 1
			elif ((c1 & 0x02) == 0):
				dx_ext0 = dx0 - 1 - SQUISH_CONSTANT_3D
				dy_ext0 = dy0 + 1 - SQUISH_CONSTANT_3D
				dz_ext0 = dz0 - 1 - SQUISH_CONSTANT_3D
				xsv_ext0 = xsb + 1
				ysv_ext0 = ysb - 1
				zsb_ext0 = zsb + 1
			else: 
				dx_ext0 = dx0 - 1 - SQUISH_CONSTANT_3D
				dy_ext0 = dy0 - 1 - SQUISH_CONSTANT_3D
				dz_ext0 = dz0 + 1 - SQUISH_CONSTANT_3D
				xsv_ext0 = xsb + 1
				ysv_ext0 = ysb + 1
				zsb_ext0 = zsb - 1

			if ((c2 & 0x01) == 0):
				dx_ext0 = dx0 + 2 - 2 * SQUISH_CONSTANT_3D
				dy_ext0 = dy0 - 2 - 2 * SQUISH_CONSTANT_3D
				dz_ext0 = dz0 - 2 - 2 * SQUISH_CONSTANT_3D
				xsv_ext0 = xsb + 2
				ysv_ext0 = ysb
				zsb_ext0 = zsb
			elif ((c1 & 0x02) == 0):
				dx_ext0 = dx0 - 2 - 2 * SQUISH_CONSTANT_3D
				dy_ext0 = dy0 + 2 - 2 * SQUISH_CONSTANT_3D
				dz_ext0 = dz0 - 2 - 2 * SQUISH_CONSTANT_3D
				xsv_ext0 = xsb
				ysv_ext0 = ysb + 2
				zsb_ext0 = zsb
			else: 
				dx_ext0 = dx0 - 2 - 2 * SQUISH_CONSTANT_3D
				dy_ext0 = dy0 - 2 - 2 * SQUISH_CONSTANT_3D
				dz_ext0 = dz0 + 2 - 2 * SQUISH_CONSTANT_3D
				xsv_ext0 = xsb
				ysv_ext0 = ysb
				zsb_ext0 = zsb + 2
			#Contribution (0,0,1)
			dx1 = dx0 - 1 - SQUISH_CONSTANT_3D
			dy1 = dy0 - 0 - SQUISH_CONSTANT_3D
			dz1 = dz0 - 0 - SQUISH_CONSTANT_3D
			attn1 = 2 - dx1 * dx1 - dy1 * dy1 - dz1 * dz1
			if (attn1 > 0): 
				attn1 *= attn1
				value = attn1 * attn1 * extrapolate(xsb + 1, ysb + 0, zsb + 0, dx1, dy1, dz1)

	 
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
	 	
	 			#Contribution (1,1,0)
			dx4 = dx0 - 1 - 2 * SQUISH_CONSTANT_3D
			dy4 = dy0 - 1 - 2 * SQUISH_CONSTANT_3D
			dz4 = dz0 - 0 - 2 * SQUISH_CONSTANT_3D
			attn4 = 2 - dx4 * dx4 - dy4 * dy4 - dz4 * dz4
			if (attn4 > 0): 
				attn4 *= attn4
				value += attn4 * attn4 * extrapolate(xsb + 1, ysb + 1, zsb + 0, dx4, dy4, dz4)
			
	 
				#Contribution (1,0,1)
			dx5 = dx4
			dy5 = dy0 - 0 - 2 * SQUISH_CONSTANT_3D
			dz5 = dz0 - 1 - 2 * SQUISH_CONSTANT_3D
			attn5 = 2 - dx5 * dx5 - dy5 * dy5 - dz5 * dz5
			if (attn5 > 0): 
				attn5 *= attn5
				value += attn5 * attn5 * extrapolate(xsb + 1, ysb + 0, zsb + 1, dx5, dy5, dz5)
				
	 
			#Contribution (0,1,1)
			dx6 = dx0 - 0 - 2 * SQUISH_CONSTANT_3D
			dy6 = dy4
			dz6 = dz5
			attn6 = 2 - dx6 * dx6 - dy6 * dy6 - dz6 * dz6
			if (attn6 > 0): 
				attn6 *= attn6
				value += attn6 * attn6 * extrapolate(xsb + 0, ysb + 1, zsb + 1, dx6, dy6, dz6)

				#First extra vertex
			attn_ext0 = 2 - dx_ext0 * dx_ext0 - dy_ext0 * dy_ext0 - dz_ext0 * dz_ext0
			if (attn_ext0 > 0):
				attn_ext0 *= attn_ext0
				value += attn_ext0 * attn_ext0 * extrapolate(xsv_ext0, ysv_ext0, zsv_ext0, dx_ext0, dy_ext0, dz_ext0)
				
		 
				#Second extra vertex
			attn_ext1 = 2 - dx_ext1 * dx_ext1 - dy_ext1 * dy_ext1 - dz_ext1 * dz_ext1
			if (attn_ext1 > 0):
				attn_ext1 *= attn_ext1
				value += attn_ext1 * attn_ext1 * extrapolate(xsv_ext1, ysv_ext1, zsv_ext1, dx_ext1, dy_ext1, dz_ext1)
	print(value)	
	return value / 18




valueArray = [[0] for x in range(128)]
i = 0

for y in range(0,128):
	for x in range(0,128):
		valueArray[y].append(eval(x, y, 0.0))
		i+=1
		print(i)
		#print(noiseValue)



