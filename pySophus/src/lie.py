'''
Created on Mar 2, 2016

@author: Francisco Dominguez
Lie group SO(2) a 2x2matrix
Lie algebra so(2) a number
'''
from visual import *
import numpy as np
import math

'''
Lie group S=(2) a 2x2matrix
'''
class SO2:
    def __init__(self,R):
        self.R=R
    def getR(self):
        return self.R
    def log(self):
        return so2.v(self.R)
    def J(self):
        return so2.G*self.R
'''
Lie algebra so(2) a twist of 1x1vector theta
'''
class so2:
    G=np.matrix([[0,-1],
                 [1, 0]])
    def __init__(self,theta):
        self.theta=theta
    def getTheta(self):
        return self.theta
    def _add_(self,op):
        R1=self.exp()
        R2=op.exp()
        R=R1*R2
        s=SO2(R)
        return s.log()
    def exp(self):
        theta=self.theta
        cs=math.cos(theta)
        sn=math.sin(theta)
        R=np.matrix([[cs,-sn],
                     [sn, cs]])
        return SO2(R)
    def x(self):
        thetax=self.theta*so2.G
        return thetax
    @staticmethod
    def v(R):
        return so2(math.atan2(R[1,0], R[0,0]))
'''
Lie group SE(2) a 3x3matrix
'''
class SE2:
    def __init__(self,T):
        self.T=T
    def T(self):
        return self.T
    def log(self):
        return se2(self.v())
    def v(self):
        w=SO2(self.T[0:2,0:2])
        t=self.T[0:2,2]
        theta=w.log().theta
        cs=math.cos(theta)
        sn=math.sin(theta)
        A=sn/theta
        B=(1-cs)/theta
        V1=1/(A**2+B**2)*np.matrix([[ A,B],
                                    [-B,A]])
        Vt=V1*t
        twist=np.matrix([[Vt[0,0]],
                         [Vt[1,0]],
                         [theta]])
        return twist
'''
Lie algebra se(2) a twist of 3x1vector vx,vy,theta
'''
class se2:
    G1=np.matrix([[0,0,1],
                  [0,0,0],
                  [0,0,0]])
    G2=np.matrix([[0,0,0],
                  [0,0,1],
                  [0,0,0]])
    G3=np.matrix([[0,-1,0],
                  [1,0,0],
                  [0,0,0]])
    def __init__(self,twist):
        self.twist=twist
    def _add_(self,op):
        R1=self.exp()
        R2=op.exp()
        R=R1*R2
        s=SE2(R)
        return s.log()
    def exp(self):
        theta=self.twist[2,0]
        cs=math.cos(theta)
        sn=math.sin(theta)
        w=so2(theta)
        t=self.twist[0:2,0]
        R=w.exp().R
        T=np.eye(3,3)
        T[0:2,0:2]=R
        V=1/theta*np.matrix([[  sn,-(1-cs)],
                             [1-cs,    sn ]])
        tmp=V*t
        T[0:2,2:3]=tmp
        return SE2(T)
    def x(self):
        return self.twist[0,0]*se2.G1+self.twist[1,0]*se2.G2+self.twist[2,0]*se2.G3
'''
Lie group SO(3) a 3x3 rotation matrix
'''
class SO3:
    def _init_(self,R):
        self.R=R
    def log(self):
        cs=(self.R.trace()-1)/2
        theta=math.acos(cs)
        sn=math.sin(theta)
        logR=theta/(2*sn)*(self.R-self.R.T)
        return so3(SO3.v(logR))
    @staticmethod
    def v(R):
        return np.matrix([[R[1,0]],
                          [R[0,2]],
                          [R[2,1]]])
    
'''
Lie algebra se(3) a 3x1matrix wx,wy,wz
'''
class so3:
    G1=np.matrix([[ 0, 0, 0],
                  [ 0, 0,-1],
                  [ 0, 1, 0]])
    G2=np.matrix([[ 0, 0, 1],
                  [ 0, 0, 0],
                  [-1, 0, 0]])
    G3=np.matrix([[ 0,-1, 0],
                  [ 1, 0, 0],
                  [ 0, 0, 0]])
    '''
    w is a 3x1matrix
    '''
    def _init_(self,w):
        self.w=w
    def _add_(self,op):
        R1=self.exp()
        R2=op.exp()
        R=R1*R2
        s=SO3(R)
        return s.log()
    def theta(self):
        return np.linalg.norm(self.w)
    def exp(self):
        wx=self.x()
        theta=np.linalg.norm(self.w)
        cs=math.cos(theta)
        sn=math.sin(theta)
        I=np.matrix([[1,0,0],
                     [0,1,0],
                     [0,0,1]])
        R=I+(sn/theta)*wx+((1-cs)/theta**2)*wx*wx
        return SO3(R)
    def x(self):
        return so3.G1*self.w[0,0]+so3.G2*self.w[1,0]+so3.G3*self.w[2,0]
    
class SE3:
    '''
    Transformation T is a 4x4matrix
    '''
    def _init_(self,T):
        self.T=T
    def log(self):
        v=np.zeros(6,1)
        R=self.T[0:2,0:2]
        w=SO3(R).v()
        
    def v(self):
        pass
class se3:
    G1=np.matrix([[ 0, 0, 0, 1],
                  [ 0, 0, 0, 0],
                  [ 0, 0, 0, 0],
                  [ 0, 0, 0, 0]])
    G2=np.matrix([[ 0, 0, 0, 0],
                  [ 0, 0, 0, 1],
                  [ 0, 0, 0, 0],
                  [ 0, 0, 0, 0]])
    G3=np.matrix([[ 0, 0, 0, 0],
                  [ 0, 0, 0, 0],
                  [ 0, 0, 0, 1],
                  [ 0, 0, 0, 0]])
    G4=np.matrix([[ 0, 0, 0, 0],
                  [ 0, 0,-1, 0],
                  [ 0, 1, 0, 0],
                  [ 0, 0, 0, 0]])
    G5=np.matrix([[ 0, 0, 1, 0],
                  [ 0, 0, 0, 0],
                  [-1, 0, 0, 0],
                  [ 0, 0, 0, 0]])
    G6=np.matrix([[ 0,-1, 0, 0],
                  [ 1, 0, 0, 0],
                  [ 0, 0, 0, 0],
                  [ 0, 0, 0, 0]])
    '''
    twist is a 6x1matrix
    '''
    def _init_(self,twist):
        self.twist=twist
    def _add_(self,op):
        R1=self.exp()
        R2=op.exp()
        R=R1*R2
        s=SE3(R)
        return s.log()
    def V(self):
        w=so3(self.twist[3:5,0])
        wx=w.x()
        theta=w.theta()
        cs=math.cos(theta)
        sn=math.sin(theta)
        I3=np.eye(3,3)
        V=I3+(1-cs)/theta**2*wx+(theta-sn)/theta**3*wx*wx
        return V
    def exp(self):
        t=self.twist[0:2,0]
        w=so3(self.twist[3:5,0])
        V=self.V()
        Vt=V*t
        T=np.eye(4,4)
        T[0:2,0:2]=w.exp()
        T[0:2,3]=Vt
        return SE3(T)
    def x(self):
        t=self.twist[0:2,0]
        w=so3(self.twist[3:5,0])
        wx=w.x()
        A=np.zeros(4,4)
        A[0:2,0:2]=wx
        A[0:2,3]=t
        return A

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    theta=-math.pi/4
    cs=math.cos(theta)
    sn=math.sin(theta)
    tx=5
    ty=5;
    T=np.matrix([[cs,-sn,tx],
                 [sn, cs,ty],
                 [0 ,  0, 1]])
    pose=SE2(T)
    print pose.T
    twist=pose.log()
    print twist
    posei=twist.exp()
    print posei.T

   