import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

def rotm2quat(R):
    quat = np.zeros((1,4));
    
    # Calculate all elements of symmetric K matrix
    K11 = R[0,0] - R[1,1] - R[2,2]
    K12 = R[0,1] + R[1,0]
    K13 = R[0,2] + R[2,0]
    K14 = R[2,1] - R[1,2]

    K22 = R[1,1] - R[0,0] - R[2,2]
    K23 = R[1,2] + R[2,1]
    K24 = R[0,2] - R[2,0]
    
    K33 = R[2,2] - R[0,0] - R[1,1]
    K34 = R[1,0] - R[0,1]

    K44 = R[0,0] + R[1,1] + R[2,2]
    
    # Construct K matrix according to paper
    K = np.array([[K11, K12, K13, K14],[K12, K22, K23, K24],
                  [K13, K23, K33, K34],[K14, K24, K34, K44]])
 
    K = K / 3
    
    eigVal, eigVec = np.linalg.eig(K)
    
    #print(eigVal)
    
    maxIdx = np.argmax(np.real(eigVal))
    
    #print(maxIdx)
    
    quat = np.real([eigVec[3,maxIdx], eigVec[0,maxIdx], eigVec[1,maxIdx], eigVec[2,maxIdx]])
    
    if quat[0] < 0:
        quat = -quat
    
    return quat


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


# Checks if a matrix is a valid rotation matrix.
def isRotationMatrix(R) :
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype = R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6

    
# Calculates rotation matrix to euler angles
def rotationMatrixToEulerAngles(R) :
 
    assert(isRotationMatrix(R))
     
    sy = math.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])
     
    singular = sy < 1e-6
 
    if  not singular :
        x = math.atan2(R[2,1] , R[2,2])
        y = math.atan2(-R[2,0], sy)
        z = math.atan2(R[1,0], R[0,0])
    else :
        x = math.atan2(-R[1,2], R[1,1])
        y = math.atan2(-R[2,0], sy)
        z = 0
        
    x = x*180/np.pi  
    y = y*180/np.pi
    z = z*180/np.pi
 
    return np.array([x, y, z])    
    
    

Data2 = sio.loadmat('C:\\Users\\REZA0\\Desktop\\Jingqiang\\Data2.mat', squeeze_me=True)

#np.save('C:\Users\REZA0\Desktop\Jingqiang\Data2.npy', Data2) 

N = 11;
V_f = 5; 
A_e =1.2;

for key, val in Data2.items():
    exec(key + '=val')
    
#x_origin = np.array([1189.8,169.5-122.5,516])
x_origin = np.array([0,0,0])

x_origin = ([1141.5 + 8.2687,140.5-94.5117,528.5+1.1513] )



x1 = np.zeros((3,len(CA_in[0][:])))   # root point
x2 = np.zeros((3,len(CA_in[0][:])))   # stub point
x3 = np.zeros((3,len(CA_in[0][:])))   #chord point

x1[0,:] = (CA_in[0][:]+x_between[0][:])/2
x1[1,:] = (CA_in[1][:]+y_between[0][:])/2
x1[2,:] = (CA_in[2][:]+intersecting[0][:])/2

x2[0,:] = CA_out[0][:]
x2[1,:] = CA_out[1][:]
x2[2,:] = CA_out[2][:]

x3[0,:] = x_between[-1][:]
x3[1,:] = y_between[-1][:]
x3[2,:] = intersecting[-1][:]


for i in range(0,len(CA_in[0][:])):
    axis_r = [0, 0, 1]
    theta_r = -np.pi/2

    [x1[0,i],x1[1,i],x1[2,i]] = np.dot(rotation_matrix(axis_r,theta_r), [x1[0,i],x1[1,i],x1[2,i]])
    [x2[0,i],x2[1,i],x2[2,i]] = np.dot(rotation_matrix(axis_r,theta_r), [x2[0,i],x2[1,i],x2[2,i]])
    [x3[0,i],x3[1,i],x3[2,i]] = np.dot(rotation_matrix(axis_r,theta_r), [x3[0,i],x3[1,i],x3[2,i]])


x1[0,:] = x1[0,:] + x_origin[0]
x1[1,:] = x1[1,:] + x_origin[1]
x1[2,:] = x1[2,:] + x_origin[2]

x2[0,:] = x2[0,:] + x_origin[0]
x2[1,:] = x2[1,:] + x_origin[1]
x2[2,:] = x2[2,:] + x_origin[2]

x3[0,:] = x3[0,:] + x_origin[0]
x3[1,:] = x3[1,:] + x_origin[1]
x3[2,:] = x3[2,:] + x_origin[2]



t=np.arange(0,len(x1[0,:]))

print(np.shape(t))

stub_len = np.zeros((len(t),))
chord_len = np.zeros((len(t),))
stub_chord_len = np.zeros((len(t),))
Angle_sc = np.zeros((len(t),))
h = np.zeros((len(t),))
x_h = np.zeros((3, len(t)))


for i in t:
    stub_len[i]=np.linalg.norm(x2[:,i]-x1[:,i])
    chord_len[i]=np.linalg.norm(x3[:,i]-x1[:,i])
    stub_chord_len[i]=np.linalg.norm(x3[:,i]-x2[:,i])
   
    Angle_sc[i] = np.arctan2(np.linalg.norm(np.cross(x2[:,i]-x1[:,i],x3[:,i]-x1[:,i])),np.dot(x2[:,i]-x1[:,i],x3[:,i]-x1[:,i]))
    
    h[i] = stub_len[i]*chord_len[i]*np.sin(Angle_sc[i])/(N*stub_chord_len[i]);
    
    x_h[:,i] = x3[:,i] + 0.5*(x2[:,i]-x3[:,i])

#plt.plot(t,Angle_sc)
#plt.show()


stub_unit_vec = np.zeros((3, len(t)))
chord_unit_vec = np.zeros((3, len(t)))

for i in t:
    stub_unit_vec[:,i]=(x2[:,i]-x1[:,i])/np.linalg.norm((x2[:,i]-x1[:,i]));
    chord_unit_vec[:,i]=(x3[:,i]-x1[:,i])/np.linalg.norm((x3[:,i]-x1[:,i]));



stub_points = np.zeros((3, N-1, len(t)))
chord_points = np.zeros((3, N-1, len(t)))
stub_chord_unit_vec = np.zeros((3, N-1, len(t)))

for i in t:
    for j in range(1,N-1):
        stub_points[:,j-1,i] = x1[:,i]+(j*stub_len[i]/N)*stub_unit_vec[:,i]
        chord_points[:,j-1,i] = x1[:,i]+(j*chord_len[i]/N)*chord_unit_vec[:,i]
        stub_chord_unit_vec[:,j-1,i]=(stub_points[:,j-1,i]-chord_points[:,j-1,i])/np.linalg.norm((stub_points[:,j-1,i]-chord_points[:,j-1,i]))

a_i = np.zeros((N-1, len(t)))
b_i = np.zeros((N-1, len(t)))
h_P = np.zeros((N-1, len(t)))
area_i = np.zeros((N-1, len(t)))

for i in t:
    for j in range(1,N-1):
        a_i[j-1,i] = (((j)/N)*stub_len[i])**2+(((j)/N)*chord_len[i])**2-2*(((j)/N)*stub_len[i])*(((j)/N)*chord_len[i])*np.cos(Angle_sc[i])
        b_i[j-1,i] = (((j+1)/N)*stub_len[i])**2+(((j+1)/N)*chord_len[i])**2-2*(((j+1)/N)*stub_len[i])*(((j+1)/N)*chord_len[i])*np.cos(Angle_sc[i])
        
        h_P[j-1,i] = stub_len[i]*chord_len[i]*np.sin(Angle_sc[i])/(N*stub_chord_len[i])/N
        
        area_i[j-1,i] = (a_i[j-1,i]+b_i[j-1,i])*h_P[j-1,i]/2


#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(x1[0,:], x1[1,:], x1[2,:], c= 'red')
#ax.scatter(x2[0,:], x2[1,:], x2[2,:], c= 'blue')
#ax.scatter(x3[0,:], x3[1,:], x3[2,:], c= 'green')
#ax.scatter(x_h[0,:], x_h[1,:], x_h[2,:], c= 'yellow')




chord_stub_subpoints = np.zeros((3, N-1, N-1, len(t)))

for i in t:
    for j in range(1,N-1):
        for k in range(1,j+1):
            if k<j+1:
                chord_stub_subpoints[:,k-1,j-1,i] = chord_points[:,j-1,i] + ((k-1)/(j+1))*(stub_points[:,j-1,i]-chord_points[:,j-1,i])
            else:
                chord_stub_subpoints[:,k-1,j-1,i] = chord_points[:,j-1,i] + ((0.5+(k-1))/(j+1))*(stub_points[:,j-1,i]-chord_points[:,j-1,i])


'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for j in range(1,N-1):
    for k in range(1,j+1):
        ax.scatter(chord_stub_subpoints[0,k,j,:],chord_stub_subpoints[1,k,j,:],chord_stub_subpoints[2,k,j,:])
'''

Plane_Normal = np.zeros((3, len(t)))
x_axis_kji = np.zeros((3, N-1, N-1, len(t)))
y_axis_kji = np.zeros((3, N-1, N-1, len(t)))
z_axis_kji = np.zeros((3, N-1, N-1, len(t)))
Center_Point_kji = np.zeros((3, N-1, N-1, len(t)))

Rotation_kji = np.zeros((3,3))
Transformation_kji = np.zeros((3, 3, N-1, N-1, len(t)))
Quaternion_kji = np.zeros((1, 4, N-1, N-1, len(t)))
Euler_kji = np.zeros((1, 3, N-1, N-1, len(t)))

Auxvar = np.zeros(3,)

for i in t:
    for j in range(1,N-1):
        for k in range(1,j+1):
            if (i>=(1/2)*len(t) and i<=(3/4)*len(t)):
                Plane_Normal[:,i]=np.cross(stub_unit_vec[:,i],chord_unit_vec[:,i])/np.linalg.norm(np.cross(stub_unit_vec[:,i],chord_unit_vec[:,i]));
                z_axis_kji[:,k-1,j-1,i]=-1.*(x_h[:,i]-x1[:,i])/np.linalg.norm(x_h[:,i]-x1[:,i])
                Center_Point_kji[:,k-1,j-1,i]=chord_stub_subpoints[:,k-1,j-1,i]
                x_axis_kji[:,k-1,j-1,i] = Plane_Normal[:,i];
                y_axis_kji[:,k-1,j-1,i]=np.cross(z_axis_kji[:,k-1,j-1,i],x_axis_kji[:,k-1,j-1,i])/np.linalg.norm(np.cross(z_axis_kji[:,k-1,j-1,i],x_axis_kji[:,k-1,j-1,i]))
                
                
                #y_axis_kji[:,k-1,j-1,i]=1.*np.cross(Plane_Normal[:,i],z_axis_kji[:,k-1,j-1,i])/np.linalg.norm(np.cross(Plane_Normal[:,i],z_axis_kji[:,k-1,j-1,i]))
                #x_axis_kji[:,k-1,j-1,i]=1.*np.cross(y_axis_kji[:,k-1,j-1,i],z_axis_kji[:,k-1,j-1,i])/np.linalg.norm(np.cross(y_axis_kji[:,k-1,j-1,i],z_axis_kji[:,k-1,j-1,i]))
                
                if np.linalg.matrix_rank(np.matrix(np.vstack((x_axis_kji[:,k-1,j-1,i],y_axis_kji[:,k-1,j-1,i],z_axis_kji[:,k-1,j-1,i]))))>2:
                    Rotation_kji = np.linalg.inv(np.matrix(np.vstack((x_axis_kji[:,k-1,j-1,i],y_axis_kji[:,k-1,j-1,i],z_axis_kji[:,k-1,j-1,i]))))
                
                Transformation_kji[:,:,k-1,j-1,i]=Rotation_kji
                
                Quaternion_kji[:,:,k-1,j-1,i]=rotm2quat(Rotation_kji)
                
                Euler_kji[:,:,k-1,j-1,i]=rotationMatrixToEulerAngles(Rotation_kji)
                
                with open('C:\\Users\\REZA0\\Desktop\\Jingqiang\\Center_Point.txt', 'ab') as output_file:
                    np.savetxt(output_file, Center_Point_kji[:,k-1,j-1,i], newline=',', fmt = '%0.4f')
                    np.savetxt(output_file, Quaternion_kji[:,:,k-1,j-1,i], newline=',', fmt = '%0.4f')
                    output_file.write(b'\r\n')


x_axis_kji = 40*x_axis_kji
y_axis_kji = 40*y_axis_kji
z_axis_kji = 40*z_axis_kji

    
'''            
            elseif (i>length(x1)/4 && i<=3*length(x1)/4)
                Plane_Normal(:,i)=cross(stub_unit_vec(:,i),chord_unit_vec(:,i))/norm(cross(stub_unit_vec(:,i),chord_unit_vec(:,i)));
                z_axis_kji(:,k,j,i)=-10.*(x_h(:,i)-x1(:,i))/norm(x_h(:,i)-x1(:,i));
                Center_Point_kji(:,k,j,i)=chord_stub_subpoints(:,k,j,i);
                y_axis_kji(:,k,j,i)=-10.*cross(Plane_Normal(:,i),z_axis_kji(:,k,j,i))/norm(cross(Plane_Normal(:,i),z_axis_kji(:,k,j,i)));
                x_axis_kji(:,k,j,i)=10.*cross(y_axis_kji(:,k,j,i),z_axis_kji(:,k,j,i))/norm(cross(y_axis_kji(:,k,j,i),z_axis_kji(:,k,j,i)));
            elseif (i>3*length(x1)/4 && i<=length(x1))
                Plane_Normal(:,i)=cross(stub_unit_vec(:,i),chord_unit_vec(:,i))/norm(cross(stub_unit_vec(:,i),chord_unit_vec(:,i)));
                z_axis_kji(:,k,j,i)=-10.*(x_h(:,i)-x1(:,i))/norm(x_h(:,i)-x1(:,i));
                Center_Point_kji(:,k,j,i)=chord_stub_subpoints(:,k,j,i);
                y_axis_kji(:,k,j,i)=10.*cross(Plane_Normal(:,i),z_axis_kji(:,k,j,i))/norm(cross(Plane_Normal(:,i),z_axis_kji(:,k,j,i)));
                x_axis_kji(:,k,j,i)=10.*cross(y_axis_kji(:,k,j,i),z_axis_kji(:,k,j,i))/norm(cross(y_axis_kji(:,k,j,i),z_axis_kji(:,k,j,i)));
'''

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x1[0,:], x1[1,:], x1[2,:], c= 'red')
ax.scatter(x2[0,:], x2[1,:], x2[2,:], c= 'blue')
ax.scatter(x3[0,:], x3[1,:], x3[2,:], c= 'green')

xo = np.array([0,0,0])
dx = np.array([30,0,0])
dy = np.array([0,30,0])
dz = np.array([0,0,30])

ax.quiver(xo[0],xo[1],xo[2],dx[0],dx[1],dx[2],color='red')
ax.quiver(xo[0],xo[1],xo[2],dy[0],dy[1],dy[2],color='green')
ax.quiver(xo[0],xo[1],xo[2],dz[0],dz[1],dz[2],color='blue')

L = 0

open('C:\\Users\\REZA0\\Desktop\\Jingqiang\\Pos_Center_Point.txt', 'w').close()
open('C:\\Users\\REZA0\\Desktop\\Jingqiang\\Rot_Center_Point.txt', 'w').close()
open('C:\\Users\\REZA0\\Desktop\\Jingqiang\\Rot_Euler_Center_Point.txt', 'w').close()

for j in range(1,2): #range(1,N-1,3):
    for i in range(np.rint((2.3/4)*len(t)).astype(int),np.rint((2.6/4)*len(t)).astype(int),5):
        for k in range(1,2):
            #print(Quaternion_kji[:,:,k-1,j-1,i])
            L = L+1 
            
            with open('C:\\Users\\REZA0\\Desktop\\Jingqiang\\Pos_Center_Point.txt', 'ab') as Pos_output_file:
                       np.savetxt(Pos_output_file, Center_Point_kji[:,k-1,j-1,i], newline=',', fmt = '%0.4f')
           
            with open('C:\\Users\\REZA0\\Desktop\\Jingqiang\\Rot_Center_Point.txt', 'ab') as Rot_output_file:
                       np.savetxt(Rot_output_file, Quaternion_kji[0,:,k-1,j-1,i], newline=',', fmt = '%0.4f')

            with open('C:\\Users\\REZA0\\Desktop\\Jingqiang\\Rot_Euler_Center_Point.txt', 'ab') as Rot_Euler_output_file:
                       np.savetxt(Rot_Euler_output_file, Euler_kji[0,:,k-1,j-1,i], newline=',', fmt = '%0.4f')
  
           
            VecStart_x = Center_Point_kji[0,k-1,j-1,i]
            VecStart_y = Center_Point_kji[1,k-1,j-1,i]
            VecStart_z = Center_Point_kji[2,k-1,j-1,i]
            
            VecEnd_x = Center_Point_kji[0,k-1,j-1,i]+x_axis_kji[0,k-1,j-1,i]
            VecEnd_y = Center_Point_kji[1,k-1,j-1,i]+x_axis_kji[1,k-1,j-1,i]
            VecEnd_z = Center_Point_kji[2,k-1,j-1,i]+x_axis_kji[2,k-1,j-1,i]
  
            ax.quiver(VecStart_x,VecStart_y,VecStart_z, VecEnd_x-VecStart_x,VecEnd_y-VecStart_y,VecEnd_z-VecStart_z,color='red')
            
            VecStart_x = Center_Point_kji[0,k-1,j-1,i]
            VecStart_y = Center_Point_kji[1,k-1,j-1,i]
            VecStart_z = Center_Point_kji[2,k-1,j-1,i]
            
            VecEnd_x = Center_Point_kji[0,k-1,j-1,i]+y_axis_kji[0,k-1,j-1,i]
            VecEnd_y = Center_Point_kji[1,k-1,j-1,i]+y_axis_kji[1,k-1,j-1,i]
            VecEnd_z = Center_Point_kji[2,k-1,j-1,i]+y_axis_kji[2,k-1,j-1,i]
  
            ax.quiver(VecStart_x,VecStart_y,VecStart_z, VecEnd_x-VecStart_x,VecStart_y-VecEnd_y,VecStart_z-VecEnd_z,color='green')
            
            VecStart_x = Center_Point_kji[0,k-1,j-1,i]
            VecStart_y = Center_Point_kji[1,k-1,j-1,i]
            VecStart_z = Center_Point_kji[2,k-1,j-1,i]
            
            VecEnd_x = Center_Point_kji[0,k-1,j-1,i]+z_axis_kji[0,k-1,j-1,i]
            VecEnd_y = Center_Point_kji[1,k-1,j-1,i]+z_axis_kji[1,k-1,j-1,i]
            VecEnd_z = Center_Point_kji[2,k-1,j-1,i]+z_axis_kji[2,k-1,j-1,i]
  
            ax.quiver(VecStart_x,VecStart_y,VecStart_z, VecEnd_x-VecStart_x,VecEnd_y-VecStart_y,VecEnd_z-VecStart_z,color='blue')

#plt.xlim(950, 1350)
#plt.ylim(-200, 200) 

plt.show()
            #print([Center_Point_kji[0,k-1,j-1,i], Center_Point_kji[0,k-1,j-1,i]+x_axis_kji[0,k-1,j-1,i]],
                    # [Center_Point_kji[1,k-1,j-1,i],Center_Point_kji[1,k-1,j-1,i]+x_axis_kji[1,k-1,j-1,i]],
                    # [Center_Point_kji[2,k-1,j-1,i],Center_Point_kji[2,k-1,j-1,i]+x_axis_kji[2,k-1,j-1,i]])

print(L)
#print(Axis)    
    
'''
soa = np.array([[0, 0, 1, 1, -2, 0], [0, 0, 2, 1, 1, 0],
                [0, 0, 3, 2, 1, 0], [0, 0, 4, 0.5, 0.7, 0]])

X, Y, Z, U, V, W = zip(*soa)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(X, Y, Z, U, V, W)
plt.show()

'''