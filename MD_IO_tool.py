import numpy as np


# In[67]:


def read_config(name):
    with open(name, 'r') as inp:
        lines = inp.readlines()
        sys_info = np.zeros((5, 2), dtype=int) # column one saves numbers, column two saves types number
        
        data = list(map(str, lines[2].split()))
        num_atom = int(data[0])
        sys_info[0,0] = num_atom 
        
        data = list(map(str, lines[3].split()))
        num_bond = int(data[0])
        sys_info[1,0] = num_bond
        
        data = list(map(str, lines[4].split()))
        sys_info[2,0] = int(data[0])
        
        data = list(map(str, lines[5].split()))
        sys_info[3,0] = int(data[0])

        data = list(map(str, lines[6].split()))
        sys_info[4,0] = int(data[0])
        
        data = list(map(str, lines[8].split()))
        sys_info[0,1] = int(data[0])
        atom_type = int (sys_info[0,1])
            
        data = list(map(str, lines[9].split()))
        sys_info[1,1] = int(data[0])
        
        data = list(map(str, lines[10].split()))
        sys_info[2,1] = int(data[0])
        
        data = list(map(str, lines[11].split()))
        sys_info[3,1] = int(data[0])
        
        data = list(map(str, lines[12].split()))
        sys_info[4,1] = int(data[0])
        
        data = list(map(str, lines[14].split()))
        boxh = float(data[1])
        
        r_pos = np.zeros((num_atom, 5))  # represents m_id, type, position
        
        for i in range(22+atom_type, 22+atom_type+num_atom):
            data = list(map(str, lines[i].split()))
            atom_id = int(data[0])-1
            r_pos[atom_id, 0] = int(data[1])
            r_pos[atom_id, 1] = int(data[2])
            for j in range(0, 3):
                r_pos[atom_id, j+2] = float(data[j+3])
                
        bond_info = np.zeros((num_bond, 3), dtype=int)
                
        for i in range(25+atom_type+num_atom, 25+atom_type+num_atom+num_bond):
            data = list(map(str, lines[i].split()))
            bond_id = int(data[0])-1    # if want to read config here, need to be -2 instead of -1, b/c there is a typo in b_id 
            for j in range(0, 3):
                bond_info[bond_id, j] = float(data[j+1])
        
        return r_pos, bond_info, sys_info, boxh


# In[68]:

def read_write(name):
    with open(name, 'r') as inp:
        lines = inp.readlines()
        sys_info = np.zeros((5, 2), dtype=int) # column one saves numbers, column two saves types number
        
        data = list(map(str, lines[2].split()))
        num_atom = int(data[0])
        sys_info[0,0] = num_atom 
        
        data = list(map(str, lines[4].split()))
        num_bond = int(data[0])
        sys_info[1,0] = num_bond
        
        sys_info[2,0] = 0 
        
        sys_info[3,0] = 0

        sys_info[4,0] = 0

        data = list(map(str, lines[3].split()))
        sys_info[0,1] = int(data[0])
        atom_type = int (sys_info[0,1])
            
        data = list(map(str, lines[5].split()))
        sys_info[1,1] = int(data[0])
        
        sys_info[2,1] = 0
        
        sys_info[3,1] = 0
        
        sys_info[4,1] = 0
        
        data = list(map(str, lines[7].split()))
        boxh = float(data[1])
        
        r_pos = np.zeros((num_atom, 5))  # represents m_id, type, position
        
        for i in range(16+atom_type, 16+atom_type+num_atom):
            data = list(map(str, lines[i].split()))
            atom_id = int(data[0]) - 1
            r_pos[atom_id, 0] = int(data[1])
            r_pos[atom_id, 1] = int(data[2])
            for j in range(0, 3):
                r_pos[atom_id, j+2] = float(data[j+3])


        vel_info = np.zeros((num_atom, 3), dtype=float)
        for i in range(19+atom_type+num_atom, 19+atom_type+num_atom+num_atom):
            data = list(map(str, lines[i].split()))
            atom_id = int(data[0]) - 1
            for j in range(0, 3):
              vel_info[atom_id, j] = data[j+1]


        bond_info = np.zeros((num_bond, 3), dtype=int)
        for i in range(22+atom_type+num_atom+num_atom, 22+atom_type+num_atom+num_atom+num_bond):
            data = list(map(str, lines[i].split()))
            bond_id = int(data[0])-1
            for j in range(0, 3):
                bond_info[bond_id, j] = float(data[j+1])
        
        return r_pos, bond_info, vel_info, sys_info, boxh


def read_rewrite_dump(dump_name, bond_name, otp_name, chl, poly_type, frame):
  
    inp = open(dump_name, 'r')
    inp_bond = open(bond_name, 'r')
    otp = open(otp_name, 'a')

    for t in range(0, frame):
      for i in range(0,3):
        line = inp.readline()

      line = inp.readline().split()
      num_atoms = line[0]

      r_pos = np.zeros((num_atoms, 5))

      line = inp.readline()

      boxh = np.zeros(3)
      for i in range(0, 3):
        line = inp.readline().split()
        boxh[i] = line[1]

      for i in range(0, num_atoms):
        line = inp.readline().split()
        atom_id = line[0]
        m_id = line[1]
        type_id = line[2]
        pos_x = line[3]
        pos_y = line[4]
        pos_z = line[5]

        r_pos[atom_id][0] = m_id
        r_pos[atom_id][1] = type_id
        r_pos[atom_id][2] = pos_x
        r_pos[atom_id][3] = pos_y
        r_pos[atom_id][4] = pos_z

      num_poly = len(np.where(r_pos[:, 1] == poly_type))
      num_bond = num_poly / chl * (chl-1)

      for i in range(0, 3):
        readline()

      line = inp_bond.readline().split()
      read_bond = line[0]

      if (read_bond != num_bond):
        print("Number of bonds are inconsistent!\n")
        exit()
      
      bond_info = np.zeros((num_bond, 2))
      for i in range(0, num_bond):
        line = inp_bond.readline().split()
        bond_info[i, 0] = line[1]
        bond_info[i, 1] = line[2]

      ## bond_shuffle



    for i in range(0, frame):
      if (i==0):
        for j in range(0, num_atoms):
          line = inp.readline().split()



# In[43]:


def bond_swap(r_pos, chl, poly_type):
    count = 0
    mole_add = 0
    atom_add = 0
    for i in range(0, np.size(r_pos[:,0])):
        if (r_pos[i,1] < poly_type):
            mole_add = max(mole_add, r_pos[i,0])  # num of NPs molecules need to added after bond swap setting
            atom_add = max(atom_add, i+1)

        elif (poly_type == r_pos[i,1]):
            atom_rand = (i-atom_add) % chl
            if (atom_rand < 0.5 * chl):
                m_id = atom_rand + 1
                if (count < 10):
                    print (m_id)
                    count += 1
            else:
                m_id = chl - atom_rand
            r_pos[i, 0] = m_id
            r_pos[i, 0] += mole_add

        else: 
            r_pos[i, 0] += 0.5 * chl
   
    print (atom_add)
    return r_pos

def read_dump(dump_name, frame):
    inp = open(dump_name, 'r')
    #n_line = np.size(inp.readlines())
    
    for t in range(0, frame):
        for i in range(0,3):
            line = inp.readline()

        line = inp.readline().split()
        num_atoms = int(line[0])
        
        if (t==0):
            r_pos = np.zeros((frame, num_atoms, 5))

        line = inp.readline()

        boxh = np.zeros(3)
        for i in range(0, 3):
            line = inp.readline().split()
            boxh[i] = line[1]
        
        line = inp.readline()

        for i in range(0, num_atoms):
            line = inp.readline().split()
            atom_id = int(line[0]) - 1
            m_id = int(line[1])
            type_id = int(line[2])
            pos_x = float(line[3])
            pos_y = float(line[4])
            pos_z = float(line[5])

            r_pos[t, atom_id, 0] = m_id
            r_pos[t, atom_id, 1] = type_id
            r_pos[t, atom_id, 2] = pos_x
            r_pos[t, atom_id, 3] = pos_y
            r_pos[t, atom_id, 4] = pos_z

    num_poly = np.size(np.where(r_pos[0, :, 1] >1))
    chl = np.size(np.where(r_pos[0, :, 0] == 2))
    #num_bond = num_poly / chl * (chl-1)
    num_chain = int(num_poly / chl)
    
    return r_pos, chl, num_chain, boxh
    

def conf_rewrite(output_file, r_pos, bond_info, sys_info, boxh):

    otp = open(output_file, "w")
    atom_num = sys_info[0, 0]
    bond_num = sys_info[1, 0]
    atom_type = sys_info[0, 1]
    
    otp.write("Generated by Entao's rewrite code\n\n")
    line = "%d atoms\n" %(atom_num)
    otp.write(line)
    line = "%d bonds\n" %(bond_num)
    otp.write(line)
    otp.write("%d angles\n" %(sys_info[2, 0]))
    otp.write('%d dihedrals\n' %(sys_info[3, 0]))
    otp.write('%d impropers\n\n' %(sys_info[4,0]))
    
    otp.write('%d atom types\n' %(atom_type))
    otp.write('%d bond types\n' %(sys_info[1, 1]))
    otp.write('%d angle types\n' %(sys_info[2, 1]))
    otp.write('%d dihedral types\n' %(sys_info[3, 1]))
    otp.write('%d improper types\n\n' %(sys_info[4, 1]))
    
    line = '%f   %f xlo xhi\n' % (-boxh, boxh)
    otp.write(line)
    line = '%f   %f ylo yhi\n' % (-boxh, boxh)
    otp.write(line)
    line = '%f   %f zlo zhi\n' % (-boxh, boxh)
    otp.write(line)
    
    otp.write('\nMasses\n\n')
    for i in range(0, atom_type):
        otp.write('%d  1.000000\n' %(i+1))
    
    otp.write('Atoms\n\n')
    for i in range(0, atom_num):
        otp.write('%d %d %d  %lf %lf %lf\n' %(i+1, int(r_pos[i,0]), int(r_pos[i,1]), r_pos[i,2], 
                                              r_pos[i, 3], r_pos[i, 4]))
    otp.write('\nBonds\n\n')
    for i in range(0, bond_num):
        otp.write('%d %d %d %d\n' %(i+1, bond_info[i, 0], bond_info[i, 1], bond_info[i, 2]))



        
'''
def reset_id(r_pos, bond_info, poly_type):
    atom_num = np.size(r_pos[:,0])
    r_pos_new = np.zeros((atom_num, 6))
    np_count = 0
    for i in range(0, atom_num):
        m_id = int(r_pos[i,0])
        a_type = int(r_pos[i,1])
        r_pos_new[i, 0] =  
        if (a_type != poly_type):
          np_count += 1
          for j in range(1, 6):
            r_pos_new[i, j] = r_pos[i, j]
    print (np_count)

    poly_count = 0
    for i in range(0, atom_num):
        m_id = int(r_pos[i,0])
        a_type = int(r_pos[i,1])
        if (a_type == poly_type):
          for j in range(0, 5):
            r_pos_new[poly_count+np_count, j] = r_pos[i, j] 
'''
