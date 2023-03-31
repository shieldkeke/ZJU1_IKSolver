from math import atan2, cos, sin, sqrt, acos, asin, pi
import numpy as np
import IK.IK as IK

class IKSolver:
    def __init__(self) -> None:
        pass
    
    def angles_to_pi(self, angles):
        for i in range(len(angles)):
            while angles[i] > pi:
                angles[i] -= 2*pi
            while angles[i] < -pi:
                angles[i] += 2*pi

    def angle_to_pi(self, angle):
        while angle > pi:
            angle -= 2*pi
        while angle < -pi:
            angle += 2*pi
        return angle
    
    def euler_angles_to_rotation_matrix(self, theta):
        R_x = [[1, 0, 0],
               [0, cos(theta[0]), -sin(theta[0])],
               [0, sin(theta[0]), cos(theta[0])]
               ]
        R_y = [[cos(theta[1]), 0, sin(theta[1])],
               [0, 1, 0],
               [-sin(theta[1]), 0, cos(theta[1])]
               ]
        R_z = [[cos(theta[2]), -sin(theta[2]), 0],
               [sin(theta[2]), cos(theta[2]), 0],
               [0, 0, 1]
               ]
        R = np.dot(R_x, np.dot(R_y, R_z))
        return R
    
    def state_to_matrix(self, state):
        R = self.euler_angles_to_rotation_matrix(state[3:])
        T = np.array([[state[0]], [state[1]], [state[2]]])
        M = np.hstack((R, T))
        M = np.vstack((M, [0, 0, 0, 1]))
        return M
    
    def solve(self, state):
        T = self.state_to_matrix(state)
        t1 = np.array([-atan2(0.023,-sqrt((T[1][3]-0.0855*T[1][2])**2 + (T[0][3]-0.0855*T[0][2])**2) - 0.023**2), 
                       -atan2(0.023,sqrt((T[1][3]-0.0855*T[1][2])**2+ (T[0][3]-0.0855*T[0][2])**2) - 0.023**2)], dtype=np.float64) + atan2(T[1][3] - 0.0855*T[1][2], T[0][3] - 0.0855*T[0][2])
        t5 = np.array([asin(T[1][2]*cos(x)-T[0][2]*sin(x)) for x in t1], dtype=np.float64)
        t5[0] = - t5[0] - np.pi
        t6 = np.array([atan2(T[1][1]*cos(x)-T[0][1]*sin(x), -T[1][0]*cos(x)+T[0][0]*sin(x)) for x in t1], dtype=np.float64)
        t6[1] = t6[1] - np.pi
        # t5 = np.array([atan2(-T[0][2]*sin(t1[i])+T[1][2]*cos(t1[i]), T[1][0]*cos(t1[i])*cos(t6[i])-T[0][0]*cos(t6[i])*sin(t1[i])-T[1][1]*cos(t1[i])*sin(t6[i])+T[0][1]*sin(t1[i])*sin(t6[i])) for i in range(2)])
        c234 = np.array([(T[0][2]*cos(t1[i])+T[1][2]*sin(t1[i]))/cos(t5[i]) for i in range(2)], dtype=np.float64)
        s234 = np.array([(T[2][0]*cos(t6[i])-T[2][1]*sin(t6[i]))/sin(t5[i]) for i in range(2)], dtype=np.float64)
        # s234 = np.array([-T[2][2]/cos(t5[i]) for i in range(2)])

        A = 0.17
        B = 0.185
        C1 = np.array([0.077*s234[i] - (T[0][3]*cos(t1[i])-0.0855 * T[0][2]*cos(t1[i]) - 0.0855 * T[1][2]*sin(t1[i]) + T[1][3]*sin(t1[i])) for i in range(2)], dtype=np.float64)
        C2 = np.array([0.077*c234[i] - (T[2][3] - 0.0855*T[2][2] - 0.23) for i in range(2)], dtype=np.float64)
        
        t2_1 = np.array([atan2(-C2[i], C1[i]) - atan2((C1[i]**2+C2[i]**2-A**2+B**2)/(2*B), sqrt(C1[i]**2+C2[i]**2-((C1[i]**2+C2[i]**2-A**2+B**2)/(2*B))**2)) for i in range(2)], dtype=np.float64)
        t2_2 = np.array([atan2(-C2[i], C1[i]) - atan2((C1[i]**2+C2[i]**2-A**2+B**2)/(2*B), -sqrt(C1[i]**2+C2[i]**2-((C1[i]**2+C2[i]**2-A**2+B**2)/(2*B))**2)) for i in range(2)], dtype=np.float64)
        t2 = np.hstack((t2_1, t2_2))
        t3_1 = np.array([atan2(-C2[i], C1[i]) - atan2((C1[i]**2+C2[i]**2+A**2-B**2)/(2*A), -sqrt(C1[i]**2+C2[i]**2-((C1[i]**2+C2[i]**2+A**2-B**2)/(2*A))**2)) - t2_1[i] for i in range(2)], dtype=np.float64)
        t3_2 = np.array([atan2(-C2[i], C1[i]) - atan2((C1[i]**2+C2[i]**2+A**2-B**2)/(2*A), sqrt(C1[i]**2+C2[i]**2-((C1[i]**2+C2[i]**2+A**2-B**2)/(2*A))**2)) - t2_2[i] for i in range(2)], dtype=np.float64)
        t3 = np.hstack((t3_1, t3_2))
        t4 = np.array([atan2(s234[i%2], c234[i%2]) - t2[i]- t3[i] for i in range(4)], dtype=np.float64)

        self.angles_to_pi(t1)
        self.angles_to_pi(t2)
        self.angles_to_pi(t3)
        self.angles_to_pi(t4)
        self.angles_to_pi(t5)
        self.angles_to_pi(t6)

        # t3_1 = np.array([acos((C1[i]**2+C2[i]**2-A**2-B**2)/(2*A*B)) for i in range(2)])
        # t3 = np.hstack((t3_1, -t3_1))
        # t2 = np.array([atan2((A*cos(t3[i])+B)*C1[i%2]-A*sin(t3[i])*C2[i%2], -(A*sin(t3[i])*C1[i%2]+(A*cos(t3[i])+B)*C2[i%2])) for i in range(4)])
        # t4 = np.array([atan2(s234[i%2], c234[i%2]) - t2[i]- t3[i] for i in range(4)])

        rst = np.array([[t1[i%2], t2[i], t3[i], t4[i], t5[i%2], t6[i%2]] for i in range(4)]).T

        return rst

if __name__ == "__main__":
    np.set_printoptions(suppress=True)
    iks = IKSolver()
    iks1 = IK.IKSolver()
    
    print("solving start state")
    print(iks1.solve([0.117, 0.334, 0.499, -2.019, -0.058, -2.190]))
    print(iks.solve([0.117, 0.334, 0.499, -2.019, -0.058, -2.190]))

    print("solving end state")
    print(iks1.solve([0.32, -0.25, 0.16, 3, 0.265, -0.84]))
    print(iks.solve([0.32, -0.25, 0.16, 3, 0.265, -0.84]))
    
    # [x, y, z, rx, ry, rz] (Cartesian, meter; X-Y'-Z'Euler, rad)
    ang_bound = np.array([[-200, 200],[-90, 90],[-120, 120],[-150, 150],[-150, 150],[-180, 180]]) /180 * np.pi
    # print(iks.euler_angles_to_rotation_matrix([1,2,3]))