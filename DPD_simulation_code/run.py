import os
import time
from multiprocessing import Pool
def setupfun(i,j,polymer,aa,cc,filename):
    if(not os.path.exists('./'+filename+'/polymer.xml')):
        os.system("mkdir " + filename);
        os.system("cp forcefield.lt ./" + filename + "/forcefield.lt")
        if (polymer):
             os.system("python3m gel_python1.py -c " + str(i) + " -o " + filename+ " -p "+ str(polymer)+ " -g "+ str(aa)+ " -u "+ str(cc));
             os.system("cp DNA_gel_test.gala ./" + filename + "/DNA_gel_test")
             os.system("cd " + filename + "&&moltemplate.sh test.lt");
             os.system("cp xml2lmpdata1.py ./" + filename + "/xml2lmpdata1.py")
             os.system("cd " + filename + "&&python2 xml2lmpdata1.py");
             os.system("cd " + filename + "&&python2 DNA_gel_test");
             os.system("cp  ./" + filename + "/test.0000050000.xml ./test.xml")
             os.system("python3m gel_python.py -c " + str(i) + " -o " + filename+ " -p "+ str(polymer)+ " -g "+ str(aa)+ " -u "+ str(cc));
             os.system("cp DNA_gel.gala ./" + filename + "/DNA_gel.gala")
        else:
            os.system("python3m gel_python.py -c " + str(i) + " -o " + filename+ " -p "+ str(polymer)+ " -g "+ str(aa)+ " -u "+ str(cc));
            os.system("cp DNA_water.gala ./" + filename + "/DNA_water.gala")
        os.system("cp gel.in ./" + filename + "/gel.in")
        os.system("sed -i \"s/0 0 5 -6 6/" + '0 0 {0} {1} {2}'.format(aa+1,-cc/2-1,cc/2+1) + "/g\" " + "./" + filename + "/gel.in");
        os.system("cp xml2lmpdata.py ./" + filename + "/xml2lmpdata.py")
        os.system("cp Rgcalcnew.py ./" + filename + "/Rgcalcnew.py")
        os.system("cp matlab.mat ./" + filename + "/matlab.mat")
        os.system("cd " + filename + "&&moltemplate.sh monomer1.lt");
        os.system("cd " + filename + "&&lmp_serial< gel.in");
        os.system("cd " + filename + "&&python2 xml2lmpdata.py polymer.data");
p = Pool(1);
polymer=1;i=5;j=4;
for aa in [2.4]:
    for cc in [2.0,2.4,2.8,3.2,3.6,4.0,4.4,4.8,5.2,5.6,6.0,6.4,6.8,7.2,7.6,8.0,8.4,8.8,9.2]:
        filename=str(aa) + "_" + str(cc)
        p.apply_async(setupfun, args=(i,j,polymer,aa,cc,filename))
p.close()
p.join()

