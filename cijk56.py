import os
import re
import math
import time
import shutil
import numpy as np
import matplotlib.pyplot as plt
import elastic_interface as elast

#### INPUTS ####
## NOTES:
# nstrain must be an even

emax    = 0.04
nstrain = 7
runjob  = 'y'
param   = 0
timer   = 2
nstep   = 1
cutoff  = 0
posfile = 'POSCAR'
poten   = 'POTCAR'
symon   = 'non'

# write generators of the point group
R4z     = np.array([    [ 0.0, 1.0, 0.0],
                        [-1.0, 0.0, 0.0],
                        [ 0.0, 0.0, 1.0]])
Rmxz    = np.array([    [ 1.0, 0.0, 0.0],
                        [ 0.0,-1.0, 0.0],
                        [ 0.0, 0.0, 1.0]])
R3111   = np.array([    [ 0.0, 1.0, 0.0],
                        [ 0.0, 0.0, 1.0],
                        [ 1.0, 0.0, 0.0]])

gen     = np.array([R4z, Rmxz, R3111])

#### SETUP ####
gamma   = np.linspace(-emax, emax , nstrain )
gdiag   = np.linspace(-emax, emax , nstrain )
h       = gdiag[1] - gdiag[0]
count   = 0
kount   = 0
C       = np.zeros((6, 6, 6))
hdir    = os.getcwd()
sig     = np.zeros((6, 6, 6, nstrain))
Cij     = np.zeros((6, 6))
coef    = elast.centraldiff(2, nstrain)
convert = 160.2177
volume  = elast.volcalc(posfile)

os.system('cp ' + posfile + ' POSCARTMP')

diff    = np.outer(coef, coef)
tol     = 1e-10
diff.real[abs(diff.real) < tol] = 0.0
ndiv    = nstrain/2
cdiff   = np.delete(diff,  (ndiv), axis=0)
Eo      = np.zeros((6,6))
dE      = np.zeros(21)
d2sde2  = np.zeros((6,21))
# load the pseuodoinverse of M
MI      = np.genfromtxt('MI.csv',delimiter=',')
M2I     = np.genfromtxt('M2I.csv',delimiter=',')
count   = 0
if param == 0:

        # determine all symmetry elements associated with the point group
        A       = elast.makegroup(gen)
        ngroup  = len(A)

        for j in xrange(6):
                for k in xrange(j,6):
                        jkdir   = 'D' + str(j) + '/' + str(k)
                        os.chdir(jkdir)

                        E       = np.zeros(nstrain)
                        for m in xrange(nstrain):
                                strdir          = str(m)
                                os.chdir(strdir)

                                F               = elast.defgrad(gamma[m], gamma[m], j+1, k+1)
                                FI              = np.linalg.inv(F)
                                J               = np.linalg.det(F)
                                E[m]            = elast.uo()
                                Cij[j,k]        = Cij[j,k] + convert/volume * E[m]*coef[m] / h**2
                                Cij[k,j]        = Cij[j,k]

                                st              = elast.getstress('OUTCAR').read()
                                s11     = st.stressx
                                s22     = st.stressy
                                s33     = st.stressz
                                s23     = st.stressyz
                                s13     = st.stressxz
                                s12     = st.stressxy
                                stress  = np.array([    [s11, s12, s13],
                                                        [s12, s22, s23],
                                                        [s13, s23, s33]])

                                s               = J * np.dot( FI, np.dot(stress, FI.T) )
                                sig[:,j,k,m]    = np.array([s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1]])

                                os.chdir('../')

                        Eo[j,k] = convert * np.dot(E,coef) / (volume*h**2)
                        plt.plot(gamma, E, '-o')
                        dE[kount]       = Eo[j,k]
                        # calculate second derivative of the stress
                        for i in xrange(6):
                                d2sde2[i,count] = np.dot(sig[i,j,k,:], coef) / h**2
                        count   = count + 1
                        kount   = kount + 1
                        os.chdir(hdir)

        # calculate the SOEC
        c2vec   = np.dot(M2I, dE)
        Cij     = np.zeros((6, 6))
        count   = 0
        for i in xrange(6):
                for j in xrange(i,6):
                        Cij[i,j]        = c2vec[count]
                        Cij[j,i]        = c2vec[count]
                        count           = count + 1

        f       = open('Cij', 'w')
        for i in xrange(6):
                f.write('{0:5f}   {1:5f}   {2:5f}   {3:5f}   {4:5f}   {5:5f}\n'.format(Cij[i,0], Cij[i,1], Cij[i,2], Cij[i,3], Cij[i,4], Cij[i,5]))

        f.close()

        print 'SOEC DONE'

        # calculate the TOEC
        nvec    = 21*6
        sigvec  = np.reshape(d2sde2.T, (nvec,))
        c3vec   = np.dot( MI, sigvec )

        f       = open('Cijk','w')
        count   = 0
        for i in xrange(6):
                for j in xrange(i,6):
                        for k in xrange(j,6):
                                f.write('c' + str(i+1) + str(j+1) + str(k+1) + ' = {0:9f}\n'.format(c3vec[count]))
                                count   = count + 1
        f.close()

        if symon[0] == 'y' or symon[0] == 'Y':
                # write the TOEC in voigt notation
                count   = 0
                cijk    = np.zeros((6,6,6))
                for i in xrange(6):
                        for j in xrange(i,6):
                                for k in xrange(j,6):
                                        cijk[i,j,k]     = c3vec[count]
                                        cijk[i,k,j]     = c3vec[count]
                                        cijk[j,i,k]     = c3vec[count]
                                        cijk[j,k,i]     = c3vec[count]
                                        cijk[k,i,j]     = c3vec[count]
                                        cijk[k,j,i]     = c3vec[count]
                                        count           = count + 1

                # symmetrize cijk into c3sym
                c3s     = np.zeros((3,3, 3,3, 3,3))
                for a in xrange(ngroup):
                  for i in xrange(3):
                    for j in xrange(3):
                      for k in xrange(3):
                        for l in xrange(3):
                          for m in xrange(3):
                            for n in xrange(3):
                              for p in xrange(3):
                                for q in xrange(3):
                                  for r in xrange(3):
                                    for s in xrange(3):
                                      for t in xrange(3):
                                        for u in xrange(3):
                                          u1    = elast.idx(p,q)
                                          u2    = elast.idx(r,s)
                                          u3    = elast.idx(t,u)
                                          c3s[i,j,k,l,m,n]= c3s[i,j,k,l,m,n] + A[a][i,p]*A[a][j,q]*A[a][k,r]*A[a][l,s]*A[a][m,t]*A[a][n,u]*cijk[u1,u2,u3]

                c3s     = c3s / ngroup
                c6s     = np.zeros((6,6,6))
                for i in xrange(3):
                 for j in xrange(3):
                  for k in xrange(3):
                   for l in xrange(3):
                    for m in xrange(3):
                     for n in xrange(3):
                      p = elast.idx(i,j)
                      q = elast.idx(k,l)
                      r = elast.idx(m,n)
                      c6s[p,q,r]= c3s[i,j,k,l,m,n]

                csymvec = np.zeros(56)
                count   = 0
                f       = open('SYMMETRIC_CIJK','w')

                for i in xrange(6):
                        for j in xrange(i,6):
                                for k in xrange(j,6):
                                        f.write('c' + str(i+1) + str(j+1) + str(k+1) + ' = {0:9f}\n'.format(c6s[i,j,k]))

                f.close()
plt.show()
#### APPLY STRAINS ####
if param == 1:
        for j in xrange(6):
                jdir    = 'D' + str(j)
                os.mkdir(jdir)
                os.chdir(jdir)
                for k in xrange(j,6):
                        kdir    = str(k)
                        os.mkdir(kdir)
                        os.chdir(kdir)

                        for m in xrange(nstrain):
                                ndir    = str(m)
                                os.mkdir(ndir)
                                os.chdir(ndir)
                                strdir  = os.getcwd()

                                # copy necessary files
                                os.chdir(hdir)
                                shutil.copy('runjob.sh' ,strdir)
                                shutil.copy('INCAR'     ,strdir)
                                shutil.copy('KPOINTS'   ,strdir)
                                shutil.copy(posfile     ,strdir)
                                shutil.copy(poten       ,strdir)
                                os.chdir(strdir)

                                # determine the deformation gradients
                                F     = elast.defgrad(gamma[m], gamma[m], j+1, k+1)
                                os.system('cp POSCAR POSCARTMP')
                                elast.writepos(F)

                                # run simulations
                                print 'mode:  {0:2d} {1:2d} {2:2d}'.format(j, k, m)

                                if runjob[0] == 'y' or runjob[0] == 'Y':
                                     # run job
                                     time.sleep(timer)
                                     os.system('sbatch runjob.sh')

                                os.chdir('../')
                        os.chdir('../')
                os.chdir('../')

