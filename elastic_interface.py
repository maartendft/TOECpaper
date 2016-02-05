import numpy as np
import os
import re
import math


def determineLength(species):
        cwd = os.getcwd()
        filename = cwd + '/' + species
        datafile = open(filename,'r')
        data = []
        for row in datafile:
                data.append(row.strip().split(','))
        return len(data)


def check(jobname):
        #make the file jobStatus
        os.system('qstat > ./tmp')
        os.system('grep ' + jobname +  ' ./tmp > jobStatus')
        os.system('rm tmp')
        n = determineLength('jobStatus')      #number of jobs remaining

        interval = 10   #number seconds before job status is checked

        while n > 0:
                import time
                time.sleep(interval)
                os.system('qstat -a > ./tmp')
                os.system('grep ' + jobname +  ' ./tmp > jobStatus')
                n = determineLength('jobStatus')
                f = open('jobsRemaining','w')
                f.write('tmp')
                f.close()
        os.system('rm jobStatus tmp jobsRemaining')

class getstress:
    def __init__(self, filename):
        self.__file__ = filename
        self.__fileopen__ = open(self.__file__)
        self.lines = self.__fileopen__.readlines()
        self.stressx = 0
        self.stressy = 0
        self.stressz = 0
        self.stressxy = 0
        self.stressyz = 0
        self.stressxz = 0
        self.finalline = ''

    def read(self):

        for line in self.lines:
            line = line.strip()
            if line.startswith("in kB"):
                self.finalline = line.split()

                self.stressx =  -0.1*float(self.finalline[2])
                self.stressy =  -0.1*float(self.finalline[3])
                self.stressz =  -0.1*float(self.finalline[4])
                self.stressxy = -0.1*float(self.finalline[5])
                self.stressyz = -0.1*float(self.finalline[6])
                self.stressxz = -0.1*float(self.finalline[7])
        return self


def writepos(F):
        ####################READ POSCAR####################
        lines = []
        poscar = "POSCARTMP"       #single input file
        f = open(poscar,'r')
        for row in f:
            lines.append((re.sub(' +',' ',row)).strip().split(' '))

        f.close()
        a0 = float(lines[1][0])
        lvstr = lines[2:5]
        if len(lines[5]) == 1:
            ntot = int(lines[5][0])
            ntype = 1
        elif len(lines[5]) == 2:
            ntot1 = int(lines[5][0])
            ntot2 = int(lines[5][1])
            ntot  = ntot1 + ntot2
            ntype = 2
        else:
            print "There are more than 2 types.  The script needs to be edited"

        posstr = lines[7:]
        #convert strings to floating points
        lv = np.zeros((3,3))
        posdir = np.zeros((ntot,3))
        pos = np.zeros((ntot,3))

        #convert the lattice vectors to floating points
        for j in xrange(3):
            for k in xrange(3):
                lv[j,k] = a0*float(lvstr[j][k])

        #convert the positions to floating points
        for j in xrange(ntot):
            for k in xrange(3):
                posdir[j,k] = float(posstr[j][k])


        ####################MAKE NEW POSCAR####################
        newlv = np.dot(lv,F)

        f = open('POSCAR','w')
        f.write('deformed\n')
        f.write('1.0000\n')

        for j in xrange(3):
            f.write('   {0:15f}   {1:15f}   {2:15f}\n'.format(newlv[j,0], newlv[j,1], newlv[j,2]))

        f.write('   ' + str(ntot) + '\n')
        f.write('Direct\n')

        for j in xrange(ntot):
            f.write('   {0:15f}   {1:15f}   {2:15f}\n'.format(posdir[j,0],posdir[j,1],posdir[j,2]))

        f.close()

def uo():
        os.system('grep "energy without" ./OUTCAR  > tmpnrg')
        os.system('tail -1 tmpnrg > tmpnrg2')
        f = open('tmpnrg2')
        tmpread = f.read()
        tmpnrg = tmpread.split()
        energy = float(tmpnrg[-1])
        os.system('rm tmpnrg tmpnrg2')
        return energy

def volcalc(poscar):
        f = open(poscar)
        tmpread = f.read()
        tmppos = tmpread.split()
        a0 = float(tmppos[1])
        v1 = np.zeros(3)
        v2 = np.zeros(3)
        v3 = np.zeros(3)
        v1[0] = a0*float(tmppos[2])
        v1[1] = a0*float(tmppos[3])
        v1[2] = a0*float(tmppos[4])
        v2[0] = a0*float(tmppos[5])
        v2[1] = a0*float(tmppos[6])
        v2[2] = a0*float(tmppos[7])
        v3[0] = a0*float(tmppos[8])
        v3[1] = a0*float(tmppos[9])
        v3[2] = a0*float(tmppos[10])
        try :
                natom = float(tmppos[11])
        except:
                natom = float(tmppos[12])
        vol =  np.dot(v1,np.cross(v2,v3))/natom
        return vol

def centraldiff(k, n):

        x       = np.linspace(-1, 1, n)*(n-1)/2
        xbar    = 0.0
        A       = np.ones((n,n))
        xrow    = x - xbar

        for i in xrange(n):
                A[i]    = xrow**(i)/np.math.factorial(i)

        b       = np.zeros(n)
        b[k]    = 1
        c       = np.linalg.solve(A , b)

        return c

def defgrad(g, h, s1, s2):
        #g:  magnitude of first strain
        #h:  magnitude of second strain
        #s1: first strain mode
        #s2: second strain mode

        F       = np.eye(3)

        ## Pure tension case
        if (s1 <=3 and s2 <=3) and s1 != s2:
                i       = s1 - 1
                j       = s2 - 1
                F[i,i]  = np.sqrt(1 + 2*g)
                F[j,j]  = np.sqrt(1 + 2*h)
        elif s1 == s2 and s1 <= 3:
                i       = s1 - 1
                j       = s2 - 1
                F[i,i]  = np.sqrt(1 + 2*g)
                F[j,j]  = np.sqrt(1 + 2*g)

        ## Pure unmixed shear case
        elif s1 == 4 and s2 == 4:
                F[1,1]  = np.sqrt(1 - 4*g**2)
                F[1,2]  = 2*g
        elif s1 == 5 and s2 == 5:
                F[0,0]  = np.sqrt(1 - 4*g**2)
                F[0,2]  = 2*g
        elif s1 == 6 and s2 == 6:
                F[0,0]  = np.sqrt(1 - 4*g**2)
                F[0,1]  = 2*g

        ## Pure unmixed shear case
        elif s1 == 4 and s2 == 4:
                F[1,1]  = np.sqrt(1 - 4*g**2)
                F[1,2]  = 2*g
        elif s1 == 5 and s2 == 5:
                F[0,0]  = np.sqrt(1 - 4*g**2)
                F[0,2]  = 2*g
        elif s1 == 6 and s2 == 6:
                F[0,0]  = np.sqrt(1 - 4*g**2)
                F[0,1]  = 2*g

        ## Pure mixed shear case
        elif (s1 == 4 and s2 == 5) or (s1 == 5 and s2 == 4):
                F[0,0]  = np.sqrt(1 - 4*g**2 - 4*h**2) / np.sqrt(1 - 4*g**2)
                F[0,1]  =-4*g*h / np.sqrt(1 - 4*g**2)
                F[0,2]  = 2*h

                F[1,1]  = (1 - 4*g**2) / np.sqrt( 1 - 4*g**2)
                F[1,2]  = 2*g

        elif (s1 == 4 and s2 == 6) or (s1 == 6 and s2 == 4):
                F[0,0]  = np.sqrt(1 - 4*g**2 - 4*h**2) / np.sqrt(1 - 4*g**2)
                F[0,1]  = 2*h / np.sqrt(1 - 4*g**2)

                F[1,1]  = np.sqrt(1 - 4*g**2)
                F[1,2]  = 2*g

        elif (s1 == 5 and s2 == 6) or (s1 == 6 and s2 == 5):
                F[0,0]  = np.sqrt(1 - 4*g**2 - 4*h**2)
                F[0,1]  = 2*h
                F[0,2]  = 2*g

        ## mixed case
        elif (s1 == 1 and s2 == 4) or (s1 == 4 and s2 == 1):
                F[0,0]  = np.sqrt(1 + 2*g)
                F[1,1]  = np.sqrt(1 - 4*h**2)
                F[1,2]  = 2*h
        elif (s1 == 1 and s2 == 5) or (s1 == 5 and s2 == 1):
                F[0,0]  = np.sqrt(1 + 2*g - 4*h**2)
                F[0,2]  = 2*h
        elif (s1 == 1 and s2 == 6) or (s1 == 6 and s2 == 1):
                F[0,0]  = np.sqrt(1 + 2*g - 4*h**2)
                F[0,1]  = 2*h

        elif (s1 == 2 and s2 == 4) or (s1 == 4 and s2 == 2):
                F[1,1]  = np.sqrt(1 + 2*g - 4*h**2)
                F[1,2]  = 2*h
        elif (s1 == 2 and s2 == 5) or (s1 == 5 and s2 == 2):
                F[1,1]  = np.sqrt(1 + 2*g)
                F[0,0]  = np.sqrt(1 - 4*h**2)
                F[0,2]  = 2*h
        elif (s1 == 2 and s2 == 6) or (s1 == 6 and s2 == 2):
                F[0,0]  = np.sqrt(1 + 2*g - 4*h**2) / np.sqrt(1 + 2*g)
                F[0,1]  = 2*h / np.sqrt(1 + 2*g)
                F[1,1]  = np.sqrt(1 + 2*g)

        elif (s1 == 3 and s2 == 4) or (s1 == 4 and s2 == 3):
                F[1,1]  = np.sqrt(1 + 2*g - 4*h**2) / np.sqrt(1 + 2*g)
                F[1,2]  = 2*h / np.sqrt(1 + 2*g)
                F[2,2]  = np.sqrt(1 + 2*g)
        elif (s1 == 3 and s2 == 5) or (s1 == 5 and s2 == 3):
                F[0,0]  = np.sqrt(1 + 2*g - 4*h**2) / np.sqrt(1 + 2*g)
                F[0,2]  = 2*h / np.sqrt(1 + 2*g)
                F[2,2]  = np.sqrt(1 + 2*g)
        elif (s1 == 3 and s2 == 6) or (s1 == 6 and s2 == 3):
                F[0,0]  = np.sqrt(1 - 4*h**2)
                F[0,1]  = 2*h
                F[2,2]  = np.sqrt(1 + 2*g)

        #strain  = 0.5*(np.dot(F, F.T) - np.eye(3))
        #tol     = 1e-10
        #strain.real[abs(strain.real) < tol] = 0.0

        return F

def belong(gcand, group):
        ng      = len(group)
        blg     = np.zeros(ng)
        for i in xrange(ng):
                matmp   = abs(gcand - group[i])
                tmp     = np.tensordot(matmp, matmp)
                if abs(tmp) < 1e-8:
                        blg[i]  = 1

        test    = int(np.sum(blg))
        if test == 0:
                return 0
        else:
                return 1

def makegroup(gen):
        ngen    = len(gen)
        # determine all the possible combinations of gen
        g1      = np.zeros((ngen**2, 3, 3))
        k       = 0
        for i in xrange(ngen):
                for j in xrange(ngen):
                        g1[k]   = np.dot(gen[i], gen[j])
                        k       = k+1

        ng1     = len(g1)
        g2      = [ g1[0] ]
        ng2     = len(g2)
        for i in xrange(1,ng1):
                # check to see if g1[i] belongs to g2
                if belong(g1[i], g2) == 0:
                        g2.append(g1[i])

                ng2     = len(g2)

        i       = 0
        ng2o    = 1
        ng2n    = 1 * ng2
        while ng2o != ng2n:
                # multiply all of the members of the group together
                g1      = 1 * g2
                ng1     = len(g1)
                k       = 0
                gtmp    = np.zeros((ng1**2,3,3))
                for j in xrange(ng1):
                        for m in xrange(ng1):
                                gtmp[k] = np.dot(g1[j], g1[m])
                                k       = k+1

                g1      = 1 * gtmp
                ng1     = len(g1)
                # get rid of duplicates
                g2      = [ g1[0] ]
                ng2     = len(g2)
                for j in xrange(1, ng1):
                        # check to see if g1[i] belongs to g2
                        if belong(g1[j], g2) == 0:
                                g2.append(g1[j])

                        ng2     = len(g2)

                ng2o    = 1 * ng2n
                ng2n    = len(g2)
                i       = i+1

        return g2

def idx(i,j):
        if i == 0 and j == 0:
                return 0
        elif i == 1 and j == 1:
                return 1
        elif i == 2 and j == 2:
                return 2
        elif (i == 1 and j == 2) or (i == 2 and j == 1):
                return 3
        elif (i == 0 and j == 2) or (i == 2 and j == 0):
                return 4
        elif (i == 0 and j == 1) or (i == 1 and j == 0):
                return 5
