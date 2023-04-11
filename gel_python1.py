import sys,getopt
import scipy.io as sio
import pymesh
print(sys.argv[1:])
opts, args = getopt.getopt(sys.argv[1:],"hc:o:p:g:u:",["chainlenarg=","outputpatharg=","polymerarg=","aa=","cc="])
#Parameters
outputpath='';
polymerarg=False
for opt, arg in opts:
  if opt == '-h':
     print('test.py -c <chainlenarg>');
     sys.exit()
  elif opt in ("-c", "--chainlenarg"):
     chainlen = int(arg)
  elif opt in ("-p", "--polymerarg"):
     polymer = int(arg)
  elif opt in ("-g", "--aa"):
     aa = float(arg)
  elif opt in ("-u", "--cc"):
     cc = float(arg)
  elif opt in ("-o", "--outputpatharg"):
     outputpath = arg


chainlen=chainlen+1;
diameterNP=30;
#aa=2;
aa=aa/2+0.25;
#cc=1;
cc=cc/2+0.25;
filename="./"+outputpath+"/test.lt"

import numpy as np
span=1.0;bn=0;Meshx=5;Meshy=5;Meshz=5;
file=open(filename, mode='w+')
file.write("import forcefield.lt\n");
file.write("\n");

b=(0.5)/(2**(1/2)/2);
NP=np.zeros([np.power(int(diameterNP*b)*4,3),3]);
cm=np.array([0,0,0]);NPnum=0;

for i in range(int(diameterNP*b)):
    for j in range(int(diameterNP*b)):
        for k in range(int(diameterNP*b)):
              NP[NPnum]=[float(i),float(j),float(k)];
              NP[NPnum+1]=[float(i)+0.5,float(j)+0.5,float(k)];
              NP[NPnum+2]=[float(i),float(j)+0.5,float(k)+0.5];
              NP[NPnum+3]=[float(i)+0.5,float(j),float(k)+0.5];
              cm=cm+[float(i),float(j),float(k)];
              cm=cm+[float(i),float(j)+0.5,float(k)+0.5];
              cm=cm+[float(i)+0.5,float(j),float(k)+0.5];
              cm=cm+[float(i)+0.5,float(j)+0.5,float(k)];
              NPnum=NPnum+4;
cm=cm/NPnum;
dott=np.subtract(NP,cm)
list=np.where( (np.sum(np.multiply(dott[:,0:2],dott[:,0:2]),axis=1) < (aa*aa)) & (np.abs(dott[:,2])<cc))
NPnew=NP[list];
NPnum=NPnew.shape[0];
cm=np.sum(NPnew,axis=0)/NPnum;
NPnew=np.subtract(NPnew,cm);
cm=np.sum(NPnew,axis=0)/NPnum;
NPnew=np.append(NPnew,cm.reshape(1,3),axis=0)
NPnew=np.append(NPnew,cm.reshape(1,3)+[0.1,0,0],axis=0)
NPnew=np.append(NPnew,cm.reshape(1,3)+[0,0.1,0],axis=0)
NPnew=np.append(NPnew,cm.reshape(1,3)+[0,0,0.1],axis=0)
NPnew=np.append(NPnew,cm.reshape(1,3)+[-0.1,0,0],axis=0)
NPnew=np.append(NPnew,cm.reshape(1,3)+[0,-0.1,0],axis=0)
NPnew=np.append(NPnew,cm.reshape(1,3)+[0,0,-0.1],axis=0)
NPnum=NPnew.shape[0];
print(NPnew)
file.write("Nanoparticle inherits PForceField {\n");
file.write("\n");
file.write("  write(\"Data Atoms\") \n");
file.write("  {\n");
for i in range(NPnum):
    file.write("    $atom:NP"+str(i)+" $mol:... @atom:NP          0.0   "+str(NPnew[i][0])+"  "+str(NPnew[i][1])+"  "+str(NPnew[i][2])+" \n");
file.write(" \n");
file.write("}\n");
file.write("}\n");


#interchain
mat = sio.loadmat('./matlab.mat');
latt = mat['lattice'];
a=(chainlen+1);
k=0;
for x in range(9):
    for y in range(9):
        for z in range(9):
            for i in range(latt.shape[0]):
                k=k+1;
                if k==1:
                    point=np.array([[(latt[i,0]+x)*a,(latt[i,1]+y)*a,(latt[i,2]+z)*a]]);
                else:
                    point=np.append(point,np.array([[(latt[i,0]+x)*a,(latt[i,1]+y)*a,(latt[i,2]+z)*a]]),axis=0);

def duplicate_removal(xy):
    if xy.shape[0] < 2:
        return xy
    _tmp = (xy*4000).astype('i4')
    _tmp = _tmp.view('S12').flatten()
    keep = np.unique(_tmp, return_index=True)[1]
    return xy[keep]
point=duplicate_removal(point)
bond=np.array([[0,0]]);
k=point.shape[0]-1;



for i in range(point.shape[0]):
    for j in range(i,point.shape[0]):
        if (np.abs(np.sum((point[i,:]-point[j,:])**2)-(chainlen+1)**2) < 1e-3):
            for p in range(chainlen):
                k=k+1;
                if (p==0):
                    point=np.append(point,[point[i,:]+(point[j,:]-point[i,:])/(chainlen+1)*(p+1)],axis=0);
                    bond=np.append(bond,np.array([[i,k]]),axis=0);
                elif(p==(chainlen-1)):
                    point=np.append(point,[point[i,:]+(point[j,:]-point[i,:])/(chainlen+1)*(p+1)],axis=0);
                    bond=np.append(bond,np.array([[k-1,k]]),axis=0);
                    bond=np.append(bond, np.array([[k, j]]), axis=0);
                else:
                    point=np.append(point,[point[i,:]+(point[j,:]-point[i,:])/(chainlen+1)*(p+1)],axis=0);
                    bond=np.append(bond,np.array([[k-1,k]]),axis=0);
minv=np.min(point,axis=0)
maxv=np.max(point,axis=0)
list1=np.unique(np.where(point==minv)[0]);
list2=np.unique(np.where(point==maxv)[0]);
listdelete=[];
for i in list2:
    m=np.where(point[i,:]==maxv)[0];
    tmp=point[i,:].copy();
    tmp[m]=minv[m];
    print(tmp)
    tmp2=np.where((point[list1,0]==tmp[0])&(point[list1,1]==tmp[1])&(point[list1,2]==tmp[2]))[0];
    print(point[list1[tmp2[0]],:],point[i,:])
    listm=(bond==i)
    bond[bond==i]=list1[tmp2[0]];
    listdelete.append(i)
bond=bond[1:,:]
print(point[np.where(abs(point[:,0]%(a*2))<1e-3)[0],0])
print(np.where(abs(point[:,0]%(a* 2)-a* 1)<1e-5)[0])
# Mesh
file.write("\n");
file.write("Mesh{\n");
for i in range(point.shape[0]):
    if(i not in listdelete):
        if((abs(point[i,0]+0.1)%(a*2)-0.1)<1e-4 and ((point[i,1]+0.1)%(a*1)-0.1)<1e-4 and ((point[i,2]+0.1)%(a*1)-0.1)<1e-4):
            file.write("     Mesh1[{0}] = new Monomer.move({1},{2},{3})\n".format(i,point[i,0],point[i,1],point[i,2]));
        elif(abs(point[i,0]%(a*2)-a*1)<1e-4 and ((point[i,1]+0.1)%(a*1)-0.1)<1e-4 and ((point[i,2]+0.1)%(a*1)-0.1)<1e-4):
            file.write("     Mesh1[{0}] = new Monomer1.move({1},{2},{3})\n".format(i,point[i,0],point[i,1],point[i,2]));
        else:
            file.write("     Mesh1[{0}] = new Monomer2.move({1},{2},{3})\n".format(i,point[i,0],point[i,1],point[i,2]));

file.write("     write(\"Data Bond List\") {\n");
for i in range(bond.shape[0]):
    bn=bn+1;
    file.write("        $bond:bb"+str(bn)+"  $atom:Mesh1[{0}]/P  $atom:Mesh1[{1}]/P\n".format(bond[i,0],bond[i,1]));
file.write("   }\n");
file.write("}\n");

#Network
file.write("write_once(\"Data Boundary\") {\n");
if(polymer):
    file.write("  {:.2f} {:.2f}  xlo xhi\n".format(-maxv[0]/2*0.47,maxv[0]/2*0.47));
    file.write("  {:.2f} {:.2f}  ylo yhi\n".format(-maxv[1]/2*0.47,maxv[1]/2*0.47));
    file.write("  {:.2f} {:.2f}  zlo zhi\n".format(-maxv[2]/2*0.47,maxv[2]/2*0.47));
else:
    file.write("  {:.2f} {:.2f}  xlo xhi\n".format(-23,23));
    file.write("  {:.2f} {:.2f}  ylo yhi\n".format(-23,23));
    file.write("  {:.2f} {:.2f}  zlo zhi\n".format(-23,23));
file.write("}\n");
file.write("NP=new Nanoparticle\n")
if(polymer):
    file.write("polymer=new Mesh.move({0},{1},{2}).scale(0.47,0.47,0.47)\n".format(-maxv[0]/2,-maxv[1]/2,-maxv[2]/2))

