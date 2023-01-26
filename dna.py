#Describle the program
print("This program reports information about DNA")
print("nucleotide sequences that may encode proteins.")

#Get the filenames
input_file = input("Input file name? ")
output_file = input("Output file name? ")

#The base bones function
#god i hope the inputs shouldnt be here, now that I think about it
def main():
    write_to_file(input_file,output_file)

def write_to_file(input_file,output_file):

    #What are the names of the sequences?
    seq_names = names(input_file)
    
    #What are the nucleotides?
    total_nucleotides = nucleotides_and_counts(input_file)
    
    nucleotides_list = total_nucleotides[0]
    nuc_counts_list = total_nucleotides[1]

    #What is the total mass?
    sums = mass(input_file)
    total_sum = sums[1]
    
    #What are the percentages?
    total_mass = percentages(input_file)

    #What is the codons list?
    total_codons = codons(input_file)
    
    #is the sequence a protein?
    total_proteins = proteins(input_file)
    

    #file1 = open(input_file, "r")
    
    #Write to the file!
    with open(output_file,"w") as file2:
            #repeat this for as many names as there are
            for line in range(len(seq_names)):
            
                file2.write(f"Region Name: {seq_names[line]}")
                file2.write(f"Nucleotides: {nucleotides_list[line]}")
                file2.write(f"Nuc. Counts: {nuc_counts_list[line]}" + "\n")
                file2.write(f"Total Mass%: {total_mass[line]} of {total_sum[line]}" + "\n")
                file2.write(f"Codons List: {total_codons[line]}" + "\n")
                file2.write(f"Is Protein?: {total_proteins[line]}" + "\n")
                file2.write("\n")

def names(input_file): #Determine list of names

    #Set up the list of names
    names = []

    #Open the file and get some names baybee!!!
    with open(input_file,"r") as file1:

        file1_contents = file1.readlines()

        #names list
        names = file1_contents[::2]
        
    #######THIS IS WHAT THIS FUNCTION SPITS OUT!!######
    return names
    
def nucleotides_and_counts(input_file): #Determine list of nucleotides
    #Set up list of nucleotides
    nucleotides = []

    #Set up list of nucleotide counts
    nuc_counts = []
    junk_count= []

    with open(input_file,"r") as file1:

        file1_contents = file1.readlines()
        
        #Determine list of nucleotides
        nucleotides = file1_contents[1::2]
        
        #turn those lil nucleotides uppercase
        for index in range(len(nucleotides)):

            nucleotides[index] = nucleotides[index].upper()
            
        #Determine nuc counts
        for sequence in nucleotides: #Do this 10 times
            
            #clear the individual nucs file
            indv_nucs_count = []
            
            #nucleotide bases/junk count
            base_A = sequence.count("A")
            base_T = sequence.count("T")
            base_G = sequence.count("G")
            base_C = sequence.count("C")
            junk = sequence.count("-")
            
            #Add the nucleotide counts
            indv_nucs_count = [base_A,base_C,base_G,base_T]

            #a third for the number of unique nucleotides (4, representing A, C, G, and T)
            nuc_counts.append(indv_nucs_count)

            #Put the junk count in the junk list
            junk_count.append(junk)

    #put \n at the end item!
    nucleotides[-1] = nucleotides[-1] + "\n"
    
    #######THIS IS WHAT THIS FUNCTION SPITS OUT!!######
    return nucleotides, nuc_counts, junk_count

def mass(input_file): #Determnes the mass of the bases
    
    #Get the nuc counts
    nucleotide_things = nucleotides_and_counts(input_file)

    #assign the return value
    nuc_counts = nucleotide_things[1]
    junk_count = nucleotide_things[2]

    #Our Lists
    masses = [135.128, 111.103, 151.128, 125.107]
    indv_nuc_masses = []
    nuc_masses = []
    sum_of_masses = []

    for count_list in nuc_counts:#Go through each list of counts

        list_sum = 0

        indv_nuc_masses = []
        
        #Loop through the numbers in the lists
        for indexs in range(4):

            #What should we multiply by to get the mass?
            current_number = count_list[indexs]
            current_mass = masses[indexs]

            #Get the mass
            nuc_mass = current_number * current_mass
        
            #put the mass in a list
            indv_nuc_masses.append(nuc_mass)

        #add them together
        list_sum = sum(indv_nuc_masses)

        #Round them
        list_sum = round(list_sum,1)
        sum_of_masses.append(list_sum)
            
        nuc_masses.append(indv_nuc_masses)

    #Add the junk and the sum_of_masses together
    for add_index in range(len(sum_of_masses)):

        correct_sum = sum_of_masses[add_index] + (junk_count[add_index]*100)
        
        sum_of_masses[add_index] = correct_sum
        
    #######THIS IS WHAT THIS FUNCTION SPITS OUT!!######
    return nuc_masses, sum_of_masses

def percentages(input_file): #Determines the percentages of the bases
    #Get the mass and sum list
    mass_and_sum = mass(input_file)

    #assign the variables
    nuc_masses = mass_and_sum[0]
    sum_of_masses = mass_and_sum[1]
    sum_index = 0

    #permanent list
    percentages = []
    
    for masses in nuc_masses:

        #start the list
        indv_percentage_list = []
        
        current_sum = sum_of_masses[sum_index]

        #Loop through the index of each list in nuc_masses. we know its ALWAYS 4 here :)
        for indexs in range(4):

            current_value = masses[indexs]

            percentage = (current_value / current_sum) *100

            percentage = round(percentage, 1)
            
            indv_percentage_list.append(percentage)
            
        percentages.append(indv_percentage_list)

        sum_index += 1

    #######THIS IS WHAT THIS FUNCTION SPITS OUT!!######
    return percentages
        
def codons(input_file): #Determines the codons for the nucleotides
    
    #get the names list
    sequence = nucleotides_and_counts(input_file)

    #assign the return value
    sequence_list = sequence[0]

    #lists
    overall_codons = []
    
    for sequence in sequence_list:
        #reset the individial codons list
        indv_codons = []

        #replaces the -'s with nothing
        fixed_sequence = sequence.replace("-","")

        #get the sequence index
        sequence_index = sequence_list.index(sequence)

        #updates the sequence list
        sequence_list[sequence_index] = fixed_sequence

        #gets the codons
        for indexs in range(0,len(fixed_sequence)-1,3):
            
            codon = fixed_sequence[indexs:indexs+3]

            #put the codons in the individial list
            indv_codons.append(codon)
            
        #put the individial codons in the overall list
        overall_codons.append(indv_codons)
        
    #######THIS IS WHAT THIS FUNCTION SPITS OUT!!######
    return overall_codons

def proteins(input_file): #Determines if the sequence is a protein
    
    #get the codons
    codons_list = codons(input_file)

    #Get the mass percentages FOR C AND G
    percentage_list = percentages(input_file)  

    #Set the variables
    start_codons = ['ATG']
    stop_codons = ['TAA', 'TAG','TGA']
    
    #min_codons = 5
    valid_length = 5

    #the percentage of mass from C and G = default of 30
    valid_percentage = 30

    #make the list
    is_protein = []

    current_index = 0
    
    #Do the calculations
    for codon in codons_list:

        current_percentage_set = percentage_list[current_index]
 
        current_index += 1

        c_percent = current_percentage_set[1]
        g_percent = current_percentage_set[2]

        #check if C and G are at least the [valid_percentage]
        if c_percent + g_percent >= valid_percentage:
            
            #check if the first codon is in the list start_codons
            if codon[0] in start_codons:
                
                #check if the last codon is in the list stop_codons
                if codon[-1] in stop_codons:

                    codons_length = len(codon)
                    
                    #check if the length of the codon is valid
                    if codons_length >= valid_length: 
                            
                        result = "Yes"

#yes I hate these else statements too but its been nearly 20 hours of work for this program
#so... it is what it is I guess
                    else:
                        result = "No"
                else:
                    result = "No"
            else:
                result = "No"
        else:
            result = "No"

        is_protein.append(result)
                
    #######THIS IS WHAT THIS FUNCTION SPITS OUT!!######   
    return is_protein

main()
