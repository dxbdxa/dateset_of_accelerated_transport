from xmlbuilder import XMLBuilder 
import numpy as np
from io import BytesIO
import re
import os
xmlfile='./polymer.xml';
lmpfile='./polymer.data';
lmp = open(lmpfile,'r+');
line=lmp.readline();
tmp=re.match("Atoms*",line);
anum=1
while(tmp==None):
    if    re.search("[0-9]+\s+atoms",line):
        natoms=int(re.search("[0-9]+",line).group())
        print(natoms)
    elif    re.search("[0-9]+\s+atom types",line):
        types=int(re.search("[0-9]+",line).group())
        print(types)
    elif    re.search("[0-9]+\s+bonds",line):
        bnum=int(re.search("[0-9]+",line).group())
        print(bnum)
    elif  re.search("[0-9]+\s+bond types",line):
        btypes=int(re.search("[0-9]+",line).group())
        print(btypes)
    elif    re.search("[0-9]+\s+angles",line):
        anum=int(re.search("[0-9]+",line).group())
        print(anum)
    elif    re.search("[0-9]+\s+angle types",line):
        atypes=int(re.search("[0-9]+",line).group())
        print(atypes)
    elif    re.search("(\-?)\d+(\.\d+)?([Ee][+-]?[\d]+)? (\-?)\d+(\.\d+)?([Ee][+-]?[\d]+)? xlo xhi",line):
        xx=re.finditer("(\-?)\d+(\.\d+)?([Ee][+-]?[\d]+)?",line)
        xlo=next(xx).group();xhi=next(xx).group();
    elif    re.search("(\-?)\d+(\.\d+)?([Ee][+-]?[\d]+)? (\-?)\d+(\.\d+)?([Ee][+-]?[\d]+)? ylo yhi",line):
        yy=re.finditer("(\-?)\d+(\.\d+)?([Ee][+-]?[\d]+)?",line)
        ylo=next(yy).group();yhi=next(yy).group();
    elif    re.search("(\-?)\d+(\.\d+)?([Ee][+-]?[\d]+)? (\-?)\d+(\.\d+)?([Ee][+-]?[\d]+)? zlo zhi",line):
        zz=re.finditer("(\-?)\d+(\.\d+)?([Ee][+-]?[\d]+)?",line)
        zlo=next(zz).group();zhi=next(zz).group();
    line=lmp.readline();
    tmp=re.search("Atoms*",line);

line=lmp.readline();
atomfull=BytesIO();
line=lmp.readline();
while(not(line=="\n")):  
    atomfull.write(line);
    line=lmp.readline();
atomfull.seek(0)
atomdata=np.loadtxt(atomfull);
atomdata=atomdata[np.argsort(atomdata[:,0])];
tmp=re.match("Bond*",line);
print(atomdata)
if(bnum):
    while(tmp==None):
        line=lmp.readline();
        tmp=re.match("Bond*",line);
    line=lmp.readline();
    bondfull=BytesIO();
    line=lmp.readline();
    while(not(line=="\n")): 
        bondfull.write(line);
        line=lmp.readline();
    bondfull.seek(0)
    bonddata=np.loadtxt(bondfull);
    bonddata=bonddata[np.argsort(bonddata[:,0])];
    

if(anum):
    try:
        tmp=re.match("Angle*",line);
        while(tmp==None):
            line=lmp.readline();
            tmp=re.match("Angle*",line);
        line=lmp.readline();
        anglefull=BytesIO();
        line=lmp.readline();
        while(not(line=="")): 
            anglefull.write(line);
            line=lmp.readline();
        anglefull.seek(0)
        angledata=np.loadtxt(anglefull);
        angledata=angledata[np.argsort(angledata[:,0])];
    except EOFError:
        pass

positiondata=atomdata[:,4:7].tolist();
PositionIO=BytesIO();
PositionIO.write("\n");
for i in range(natoms):
    PositionIO.write("   {0}   {1}   {2}\n".format(positiondata[i][0],positiondata[i][1],positiondata[i][2]));
PositionIO.seek(0)
text1=PositionIO.read();
PositionIO.close();

print(np.shape(atomdata)[1])
if(np.shape(atomdata)[1]>7):
    imagedata=atomdata[:,7:10].tolist();
    ImageIO=BytesIO();
    ImageIO.write("\n");
    for i in range(natoms):
        ImageIO.write("{0} {1} {2}\n".format(imagedata[i][0],imagedata[i][1],imagedata[i][2]));
    ImageIO.seek(0)
    text2=ImageIO.read();
    ImageIO.close();

MassIO=BytesIO();
MassIO.write("\n");
for i in range(natoms):
    MassIO.write("{0}\n".format("1.00"));
MassIO.seek(0)
text3=MassIO.read();
MassIO.close();

typedata=atomdata[:,2]
TypeIO=BytesIO();
TypeIO.write("\n");
for i in range(natoms):
    TypeIO.write("{0}\n".format(chr(int(typedata[i])+64)));
TypeIO.seek(0)
text4=TypeIO.read();
TypeIO.close();
if(bnum):
    BondIO=BytesIO();
    bondtype,blo,blm=np.unique(bonddata[:,1],return_index=True,return_inverse=True);
    def findtype(a,b):
        return int(atomdata[int(bonddata[blo[int(a)-1],b])-1,2])
    bondtype=bondtype.tolist();
    bondmap=map(lambda x: [findtype(bondtype.index(x),2),findtype(bondtype.index(x),3)],bondtype);
    print(bondmap)
    BondIO.write("\n");
    for i in range(bnum):
        BondIO.write("{0}-{1} {2} {3}\n".format(chr(int(atomdata[int(bonddata[i,2])-1,2])+64),chr(int(atomdata[int(bonddata[i,3])-1,2])+64),int(bonddata[i,2]-1),int(bonddata[i,3]-1)));
        if(i<5):
            print(atomdata[int(bonddata[i,2])-1,2],atomdata[int(bonddata[i,3])-1,2])
    BondIO.seek(0)
    text5=BondIO.read();
    BondIO.close();


if(anum):
    AngleIO=BytesIO();
    AngleIO.write("\n");
    for i in range(anum):
        AngleIO.write("{0}-{1}-{2} {3} {4} {5}\n".format(chr(int(atomdata[int(angledata[i,2])-1,2])+64),chr(int(atomdata[int(angledata[i,3])-1,2])+64),chr(int(atomdata[int(angledata[i,4])-1,2])+64),int(angledata[i,2]-1),int(angledata[i,3]-1),int(angledata[i,4]-1)));
    AngleIO.seek(0)
    text6=AngleIO.read();
    AngleIO.close();
    
body=1;
    
if(body):
    BodyIO=BytesIO();
    BodyIO.write("\n");
    def findbodytype(a):
        if(atomdata[i,2]==4):
            return 0
        else:
            return -1
    for i in range(natoms):
        BodyIO.write("{0}\n".format(findbodytype(i)));
    BodyIO.seek(0)
    text7=BodyIO.read();
    BodyIO.close();
orientation=1;
if(orientation):
    orientationIO=BytesIO();
    orientationIO.write("\n");
    def findbodytype(a):
        if(atomdata[i,2]==4):
            return 1
        else:
            return 0
    for i in range(natoms):
        orientationIO.write("0.00 0.00 {0}\n".format(findbodytype(i)));
    orientationIO.seek(0)
    text8=orientationIO.read();
    orientationIO.close();

npcm=np.where((atomdata[typedata==4,4:7]**2).sum(axis=1)<1e-5)[0];
npnumber=np.where(typedata==4)[0];
npatoms=npnumber.shape[0];
print(npatoms)
npcm=npatoms-7;
vv=npatoms-7;
xv=npatoms-6;
yv=npatoms-5;
zv=npatoms-4;

np.savetxt('npnumber.txt',np.array([npcm,vv,xv,yv,zv]).astype('int'),fmt='%d')

x = XMLBuilder('galamost_xml',version='1.3') 
y=x.configuration(natoms=str(natoms),time_step='0',dimensions="3")
y.box(lx=str(float(xhi)-float(xlo)),ly=str(float(yhi)-float(ylo)),lz=str(float(zhi)-float(zlo)))
y.position(text1,num=str(natoms));
if(np.shape(atomdata)[1]>7):
    y.image(text2,num=str(natoms));
y.mass(text3,num=str(natoms));
y.type(text4,num=str(natoms));
y.body(text7,num=str(natoms));
y.orientation(text8,num=str(natoms));
if(bnum):
    y.bond(text5,num=str(bnum));
if(anum):
    y.angle(text6,num=str(anum));

etree_node = ~x # <= return xml.etree.ElementTree object 
xml=open(xmlfile,'w+')
xml.write(str(x)); 
xml.write("\n")
xml.close();
