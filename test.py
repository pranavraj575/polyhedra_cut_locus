import numpy as np

norm = np.linalg.norm

# Creating points
p1 = np.array([0,-4/3])
p2 = np.array([2, 0])
p3 = np.array([5, 6])
p2=np.expand_dims(p2,0)
p3=np.expand_dims(p3,0)
p1=np.expand_dims(p1,0)
p2=np.expand_dims(p2,0)
p3=np.expand_dims(p3,0)
# Display original points
print("p1:\n",p1,"\n")
print("p2:\n",p2,"\n")
print("p3:\n",p3,"\n")
print(np.cross(p2-p1, p1-p3,axis=-1))
# Finding distance of p3 from a
# line connecting p1 and p2
print(norm(np.cross(p2-p1, p1-p3,axis=-1),axis=-1))
res = norm(np.cross(p2-p1, p1-p3,axis=-1),axis=-1)/norm(p2-p1,axis=-1)
print(np.array([True,True]) + np.array([True,False]))
# Display result
print("Result:\n",res,"\n")