
#Written by Jonas Kolker and Andres Cardenas
#December 2017

#This program models protein folding on a basic, 2D level
#It can both create proteins from scratch as well as go off of premade data files

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as anim
import argparse
import time

#########################
######## CLASSES ########
#########################

#takes in a base point and length - has a second point attribute that should be assigned to an Atom
class Bond(object):

	def __init__(self, basePoint, length):
		self.basePoint = basePoint
		self.secPoint  = np.array([(basePoint[0] + length), basePoint[1]])
		self.length    = length
		
	def radianTake(self, radian):
		legX = self.length * np.cos(radian)
		legY = self.length * np.sin(radian)
		self.secPoint = np.array([legX + self.basePoint[0], legY + self.basePoint[1]])

#takes in a charge and position as arqument
class Atom(object):
	def __init__(self, charge, position):
		self.pos  = position
		self.q    = charge

#########################
###### Functions ########
#########################

def distanceGet(point1, point2):
	squareRootMe = (point1[0] - point2[0])**2 + (point1[1] - point2[1])**2
	distance     = np.sqrt(squareRootMe)
	return distance

def rotate(vector, radians):
	rotation 	= np.array([ [np.cos(radians), - np.sin(radians)],
						     [np.sin(radians),   np.cos(radians)] ])
	out 		= rotation.dot(vector)
	return out

def pairEnrg(particle1, particle2):
	energy  = (alpha * particle1.q * particle2.q) / distanceGet(particle1.pos, particle2.pos)
	return energy

#particle must be an array of Atoms
def totalEnergyGet(particle, matrix):
	for i in range(len(particle)):
		for j in range(i+1, len(particle)):
			if np.abs(i-j) != 1:
				matrix[i,j] = pairEnrg(particle[i], particle[j])
			else:
				matrix[i,j] = bond_Enrg
	totalEnrg = (np.sum(matrix))#deleted the divided by two bc it's unnecesary - see the i+1
	return totalEnrg

#creates the bonds at angles and places Atoms at each end
def bondParticleCreator(particle, bonds, numParticles):
	for i in range(numParticles-1):
		#random radian between range between pi/4 and -pi/4
		randRad = np.pi/4
		r = particle[i].pos
		if i%2:
			#creates a bond with a base point at the previous Atom with randRad radian
			bonds.append(Bond(r, bond_length))
			bonds[i].radianTake(randRad)
			#sets the position of the next Atom as the second point of the just-made bond
		else:
			#does same as above but with a negative randRad
			bonds.append(Bond(r, bond_length))
			bonds[i].radianTake(-randRad)
		particle[i+1].pos = bonds[i].secPoint

def ranProtTurn(particlePos_temp, randPoint, randRad, i):
	tempPos = rotate(particlePos_temp[i].pos - randPoint, randRad)
	particlePos_temp[i].pos = tempPos + randPoint
	return

def fileWrite(fname, array, particle):
	for i in range(numParticles):
		array[i,0] = particle[i].q
		array[i,1] = particle[i].pos[0]
		array[i,2] = particle[i].pos[1]
	np.save(fname, array)
	return

def fileRead(fname, numParticles):
	particle = []
	inputArray = np.load(fname + '.npy')
	
	for i in range(numParticles):
		#print(i)
		particle.append(Atom(inputArray[(i,0)], (inputArray[(i,1)], inputArray[(i,2)]) ))
	return particle

def constantGet(fname):
	constants 		= open(fname, 'r')
	vals 			= []

	for line in constants:
		 vals.append(float(line.split()[-1]))

	numParticles 	= int(vals[0])
	bond_length		= vals[1]
	num_iteration	= int(vals[2])
	alpha			= vals[3]
	T				= vals[4]
	elecCharge		= vals[5]
	k 				= vals[6]
	bond_Enrg		= vals[7]


	return 	numParticles, bond_length, num_iteration, alpha, T, elecCharge, k, bond_Enrg

def tempCheck(temp, energy0, energy1):
	R       = np.random.random()

	#this is e to the negative delta of energy divided by kT. This will be a value between 0 and 1: a probability. 
	#this will then be compared to the random probability number that is generated. The higher the temperature, the more
	#likely that checkMe will be greater than the constant R
	
	checkMe = np.exp(-(energy0 - energy1)/kT)

	
	print('\ncheckMe = {}'.format(checkMe))
	print('kT = {}'.format(kT))
	print('energyDif = {} \n'.format(energy0 - energy1))

	if min(1., checkMe) >= R:
		print('Temp was accounted for')
		return True
	elif min(1., checkMe) < R:

		return False

def makeBaseProtein(particle, inputMode):
	if inputMode == False:
		for i in range(numParticles):
			if i % 2 == 0:
				particle.append(Atom(elecCharge, np.array([0.,0.]))) 
				print(particle[i].q)
			else:
				particle.append(Atom(-elecCharge, np.array([0.,0.])))
				print(particle[i].q)
		
		bonds 	     = []
		bondParticleCreator(particle, bonds, numParticles)
	
	elif inputMode:
		particle = fileRead(fname, numParticles)
		print('input mode')	

	return particle, bonds


#########################
##### MAIN DEF ##########
#########################

def main():

	#---------------------------------------------------------------------------------
	#Get constants from the np file and converts them to floats
	
	constants 		= constantGet('./constants.par')
	
	global numParticles
	global bond_length
	global num_iteration
	global alpha
	global T
	global elecCharge
	global k
	global bond_Enrg
	global kT

	numParticles 	= constants[0]
	bond_length		= constants[1]
	num_iteration 	= constants[2]
	alpha			= constants[3]
	T 				= constants[4]
	elecCharge		= constants[5]
	k 				= constants[6]
	bond_Enrg		= constants[7]
	kT 				= k * T

	#create a binary file to store the protein
	
	global fname
	fname 			= './dataOutput'
	data 			= np.zeros((numParticles, 3))

	#create switches - parameters in the bash execution of the code
	
	parser = argparse.ArgumentParser(description= 'This is a 2D simple protein folding code.')
	parser.add_argument('-i', action='store_true', dest='inputMode', help='Allows the code to read a pre-existing protein file')
	inputs = parser.parse_args()
	inputMode = inputs.inputMode

	#---------------------------------------------------------------------------------
	##################################################################################
	#---------------------------------------------------------------------------------
	
	#Creates the arrays of particles (with different charges), particle positions, and bonds
	
	particle, bonds 		= [], []
	particle, bonds 	 	= makeBaseProtein(particle, inputMode)
	
	#---------------------------------------------------------------------------------
	##################################################################################
	#---------------------------------------------------------------------------------
	
	#Creates a canvas, creates an axis on a 1 supblot in the first row and first column (no autoscaling)
	#Sets the x and y axis limits to the number of bonds divided by the bond lengthcd ..

	fig = plt.figure()
	ax = fig.add_subplot(111, autoscale_on = False)
	ax.set_xlim(- (numParticles - 1 ) * (.5 * bond_length), (numParticles - 1 ) * (.5 * bond_length))
	ax.set_ylim(- (numParticles - 1 ) * bond_length, (numParticles - 1 ) * bond_length)
	
	#---------------------------------------------------------------------------------
	##################################################################################
	#---------------------------------------------------------------------------------
	
	#For each iteratoin the any point after the choosen point will be rotated by a random radian 
	#Also finds the energy between two points 
	
	energyMatrix 	= np.zeros([numParticles, numParticles])
	pastTotalEnergy = totalEnergyGet(particle, energyMatrix)

	print('Starting energy: ', pastTotalEnergy)
	
	#---------------------------------------------------------------------------------
	##################################################################################
	#---------------------------------------------------------------------------------	
	
	#for each step, all of the points after the choosen particle are rotated and their 
	#energies are compared to the previous
	acceptedSteps = 0
	for step in range(num_iteration):	
		particleTemp= particle.copy()
		
		randIdx    = np.random.random_integers(0, numParticles -1)
		randPoint  = particleTemp[randIdx].pos
		randRad    = (np.pi/64) * np.random.random_sample() - (np.pi/128)

		for i in range(randIdx+1, len(particle)):
			#creates a temporary position that is equivalent to the current point minus the vertex of the choosen point
			#rotates the temporaray position by the random radian
			#sets the current point to the temporary point
			ranProtTurn(particleTemp, randPoint, randRad, i)
		# Accept this step if energy is explicitly less and copy the 
		#temporary positions into the current positions
		totalEnergy = totalEnergyGet(particleTemp, energyMatrix)
		if acceptedSteps % 200 == 0:
				fileWrite(fname, data, particle)
				print(step, acceptedSteps, totalEnergy)
		print('pastTotalEnergy: {}'.format(pastTotalEnergy))
		print('totalEnergy: {}'.format(totalEnergy))

		# Accept this step if the energy is less w/ temp and 
		#copy the temporary positions into the current positions
		if tempCheck(T, pastTotalEnergy, totalEnergy):
			pastTotalEnergy = totalEnergy
			particle = particleTemp.copy()
			acceptedSteps += 1
		print('pastTotalEnergy: {}'.format(pastTotalEnergy))
		print('totalEnergy: {}'.format(totalEnergy))
	
	#---------------------------------------------------------------------------------
	##################################################################################
	#---------------------------------------------------------------------------------
	
		#creates a numpy array of the particles' positions with x and y values assigned accordingly
		particlePos = np.array([particle[i].pos for i in range(numParticles)])

		x, y = particlePos.T
		ax.clear()
		ax.plot(x,y)
		ax.set_xlim(- (numParticles - 1 ) * bond_length, (numParticles - 1 ) * bond_length)
		ax.set_ylim(- 20 * bond_length, 20 * bond_length)
		ax.set_title('Step: {}, E = {}' .format(step, pastTotalEnergy))
		plt.pause(0.1)



	plt.show()


#########################
##### CALLING MAIN ######
#########################

if __name__ == "__main__":
	main()
