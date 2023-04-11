import sys,getopt
import scipy.io as sio
import pymesh
import  xml.dom.minidom
from io import StringIO
import networkx as nx
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

npfile='{0}_{1}'.format(aa,cc)
chainlen=chainlen+1;
diameterNP=30;sca=0.47;
#aa=2;
aa=aa/2;
#cc=1;
cc=cc/2;
filename="./"+outputpath+"/monomer1.lt"

import numpy as np
span=1.0;bn=0;Meshx=5;Meshy=5;Meshz=5;
file=open(filename, mode='w+')
file.write("import forcefield.lt\n");
file.write("\n");

#b=(0.5)/(2**(1/2)/2);
#NP=np.zeros([np.power(int(diameterNP*b)*4,3),3]);
cm=np.array([0,0,0]);NPnum=0;

#Nanoparticle
# for i in range(int(diameterNP*b)):
    # for j in range(int(diameterNP*b)):
        # for k in range(int(diameterNP*b)):
              # NP[NPnum]=[float(i),float(j),float(k)];
              # NP[NPnum+1]=[float(i)+0.5,float(j)+0.5,float(k)];
              # NP[NPnum+2]=[float(i),float(j)+0.5,float(k)+0.5];
              # NP[NPnum+3]=[float(i)+0.5,float(j),float(k)+0.5];
              # cm=cm+[float(i),float(j),float(k)];
              # cm=cm+[float(i),float(j)+0.5,float(k)+0.5];
              # cm=cm+[float(i)+0.5,float(j),float(k)+0.5];
              # cm=cm+[float(i)+0.5,float(j)+0.5,float(k)];
              # NPnum=NPnum+4;
#mesh=pymesh.generate_cylinder([0,0,-cc], [0,0,cc], aa, aa, num_segments=16)
#mesh=pymesh.generate_icosphere(aa,[0,0,0],1)
#mesh=pymesh.split_long_edges(mesh, 0.4)[0]
#mesh=pymesh.collapse_short_edges(mesh, abs_threshold=0.05, rel_threshold=True, preserve_feature=False)[0]
mat = sio.loadmat('./'+npfile+'.mat');
#NP=mesh.vertices;
NP=mat['p'];
NPnum=NP.shape[0];
cm=cm/NPnum;
dott=np.subtract(NP,cm)
NPnew=NP;
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
k=0;p=0;
for x in range(9):
    for y in range(9):
        for z in range(9):
            for i in range(latt.shape[0]):
                k=k+1;
                if k==1:
                    point=np.array([[(latt[i,0]+x)*a,(latt[i,1]+y)*a,(latt[i,2]+z)*a]]);
                else:
                    point=np.append(point,np.array([[(latt[i,0]+x)*a,(latt[i,1]+y)*a,(latt[i,2]+z)*a]]),axis=0);
#gap=np.zeros((0,6));
#for x in range(6):
    #for y in range(6):
        #for z in range(6):
            #gap=np.append(gap,np.array([x,y,z])*a+np.array([[0,0,0.5],[0,0.5,0],[0.5,0,0],[0.5,0.5,0.5],[0.25,0.25,0.75],[0.25,0.75,0.25],[0.75,0.25,0.25],[0.75,0.75,0.75]])*a,axis=0)

def duplicate_removal(xy):
    if xy.shape[0] < 2:
        return xy
    _tmp = (xy*4000).astype('i4')
    _tmp = _tmp.view('S12').flatten()
    keep = np.unique(_tmp, return_index=True)[1]
    return xy[keep]
point=duplicate_removal(point)
minv=np.min(point,axis=0)
maxv=np.max(point,axis=0)
#r=np.zeros((gap.shape[0],point.shape[0],3))
#for j in range(point.shape[0]):
    #for i in range(gap.shape[0]):
        #r[i,j,:]=(point[j,:]-gap[i,:]);
#r[np.where(abs(r)>(maxv[0]/2))]=maxv[0]-abs(r[np.where(abs(r)>(maxv[0]/2))])
#r=np.linalg.norm(r, axis=2);
#gappoint=np.where(abs(r-(a*0.5))<1e-3)
#gappair=np.array([gappoint[0],gappoint[1]]).T
bond=np.array([[0,0]]);chains=np.zeros((0,chainlen+2)).astype('int');
k=point.shape[0]-1;
for i in range(point.shape[0]):
    for j in range(i,point.shape[0]):
        if (np.abs(np.sum((point[i,:]-point[j,:])**2)-(chainlen+1)**2) < 1e-3):
            chainstmp=[];
            for p in range(chainlen):
                k=k+1;
                if (p==0):
                    point=np.append(point,[point[i,:]+(point[j,:]-point[i,:])/(chainlen+1)*(p+1)],axis=0);
                    bond=np.append(bond,np.array([[i,k]]),axis=0);
                    chainstmp.append(i)
                    chainstmp.append(k)
                elif(p==(chainlen-1)):
                    point=np.append(point,[point[i,:]+(point[j,:]-point[i,:])/(chainlen+1)*(p+1)],axis=0);
                    bond=np.append(bond,np.array([[k-1,k]]),axis=0);
                    bond=np.append(bond, np.array([[k, j]]), axis=0);
                    chainstmp.append(k)
                    chainstmp.append(j)
                else:
                    point=np.append(point,[point[i,:]+(point[j,:]-point[i,:])/(chainlen+1)*(p+1)],axis=0);
                    bond=np.append(bond,np.array([[k-1,k]]),axis=0);
                    chainstmp.append(k)
            chains=np.append(chains, np.array([chainstmp]), axis=0);
minv=np.min(point,axis=0)
maxv=np.max(point,axis=0)

list1=np.unique(np.where(point==minv)[0]);
list2=np.unique(np.where(point==maxv)[0]);
listdelete=[];
for i in list2:
    m=np.where(point[i,:]==maxv)[0];
    tmp=point[i,:].copy();
    tmp[m]=minv[m];
    tmp2=np.where((point[list1,0]==tmp[0])&(point[list1,1]==tmp[1])&(point[list1,2]==tmp[2]))[0];
    bond[bond==i]=list1[tmp2[0]];
    chains[chains==i]=list1[tmp2[0]]
    #gappair=gappair[(gappair[:,1]!=i)]
    listdelete.append(i)
bond=bond[1:,:]
#gaplist=gappair[:,1].reshape(3*3*3*8,6)
#cages=np.empty((0,0));k=0;
#for i in range(gaplist.shape[0]):
    #tmp=np.empty(0).astype('int');
    #for j in range(gaplist.shape[1]):
        #tmp1=np.where(chains[:,[0,chainlen+1]]==gaplist[i,j])[0];
        #tmp2=chains[tmp1][:,[0,chainlen+1]];
        #tmp3=point[tmp2]-gap[i]
        #tmp3[np.where(abs(tmp3)>(maxv[0]/2))]=maxv[0]-abs(tmp3[np.where(abs(tmp3)>(maxv[0]/2))])
        #tmp1=np.delete(tmp1,np.where(abs(tmp3)>(a*0.6))[0])
        #tmp=np.append(tmp,tmp1);
    #if(k==0):
        #cages=np.array([np.unique(chains[tmp])]);
    #else:
        #cages=np.append(cages,np.array([np.unique(chains[tmp])]),axis=0)
    #k=k+1;
    
#print(r[57])
xmlfile="./"+outputpath+"/test.0000050000.xml";
dom = xml.dom.minidom.parse(xmlfile)
def get_nodevalue(node, index = 0):
    return node.childNodes[index].nodeValue if node else ''

root = dom.documentElement
configurtion=time_step=root.getElementsByTagName("configuration");
time_step=root.getElementsByTagName("configuration")[0].attributes.getNamedItem("time_step").nodeValue
natoms=root.getElementsByTagName("configuration")[0].attributes.getNamedItem("natoms").nodeValue

position=configurtion[0].getElementsByTagName("position")

positionfile=StringIO();
positionfile.write(position[0].firstChild.data)
positionfile.seek(0)
pos=np.loadtxt(positionfile)

# Mesh
file.write("\n");
file.write("Mesh{\n");
k=0;
for i in range(point.shape[0]):
    if(i not in listdelete):
        
        if((abs(point[i,0]+0.1)%(a*2)-0.1)<1e-4 and ((point[i,1]+0.1)%(a*1)-0.1)<1e-4 and ((point[i,2]+0.1)%(a*1)-0.1)<1e-4):
            file.write("     Mesh1[{0}] = new Monomer.move({1},{2},{3})\n".format(i,pos[k,0],pos[k,1],pos[k,2]));
        elif(abs(point[i,0]%(a*2)-a*1)<1e-4 and ((point[i,1]+0.1)%(a*1)-0.1)<1e-4 and ((point[i,2]+0.1)%(a*1)-0.1)<1e-4):
            file.write("     Mesh1[{0}] = new Monomer1.move({1},{2},{3})\n".format(i,pos[k,0],pos[k,1],pos[k,2]));
        else:
            file.write("     Mesh1[{0}] = new Monomer2.move({1},{2},{3})\n".format(i,pos[k,0],pos[k,1],pos[k,2]));
        #cages[cages==i]=k;
        #gaplist[gaplist==i]=k;
        k=k+1;
#cagedata={'cages': cages,'gaplist':gaplist}
#np.save("./"+outputpath+"/cage.npy",cagedata)

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
    file.write("polymer=new Mesh\n")

