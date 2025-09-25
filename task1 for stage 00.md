
# **BASh Basic & Bioinformatics Software Installation**

This repository contains exercises and commands related to basic shell operations (Bash) and bioinformatics software installation.

---

## **Project 1: BASh Basic**

### **1. Print Your Name**
Print your name to the terminal using the following command:

```bash
habiba_faried22@cloudshell:~$ echo "habiba"
```

### **2. Create a Folder Titled Your Name**
Create a directory with your name.

```bash
habiba_faried22@cloudshell:~$ mkdir habiba
```

### **3. Create a Directory Titled `biocomputing` and Change to That Directory**
Use a single line of command to create the `biocomputing` directory and navigate into it.

```bash
habiba_faried22@cloudshell:~$ mkdir biocomputing && cd biocomputing
```

### **4. Download Files**
Download three required files using `wget`.

```bash
habiba_faried22@cloudshell:~/biocomputing$ wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
habiba_faried22@cloudshell:~/biocomputing$ wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
habiba_faried22@cloudshell:~/biocomputing$ wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
```

### **5. Move the `.fna` File to the Folder Titled Your Name**
Move the `wildtype.fna` file to the folder created in step 2.

```bash
habiba_faried22@cloudshell:~/biocomputing$ mv wildtype.fna ../habiba/
```

### **6. Delete the Duplicate `.gbk` File**
Delete the duplicate `wildtype.gbk` file.

```bash
habiba_faried22@cloudshell:~/biocomputing$ rm wildtype.gbk.1
```

### **7. Confirm if the `.fna` File is Mutant or Wild Type**
Search for `tatatata` in the `.fna` file to determine if it is a mutant.

```bash
habiba_faried22@cloudshell:~/habiba$ if grep -q "tatatata" ../habiba/wildtype.fna; then
    echo "Mutant"
elif grep -q "tata" ../habiba/wildtype.fna; then
    echo "Wild type"
else
    echo "Pattern not found"
fi
```

### **8. If Mutant, Print All Matching Lines into a New File**
If it’s a mutant, save all matching lines to a file called `mutant_lines.txt`.

```bash
habiba_faried22@cloudshell:~/habiba$ grep 'tatatata' ../habiba/wildtype.fna > mutant_lines.txt
```

### **9. Count the Number of Lines in the `.gbk` File**
Count the total number of lines (excluding the header) in the `.gbk` file.

```bash
habiba_faried22@cloudshell:~/biocomputing$ wc -l wildtype.gbk
```

### **10. Print the Sequence Length of the `.gbk` File**
Print the sequence length of the `.gbk` file using the `LOCUS` tag.

```bash
habiba_faried22@cloudshell:~/biocomputing$ awk 'NR==1{for(i=1;i<=NF;i++) if($(i+1)=="bp" && $i ~ /^[0-9]+$/){print $i; exit}}' wildtype.gbk
```

### **11. Print the Source Organism of the `.gbk` File**
Print the **SOURCE** tag from the `.gbk` file.

```bash
habiba_faried22@cloudshell:~/biocomputing$ awk '/^SOURCE/{sub(/^SOURCE[ 	]+/, ""); print; exit}' wildtype.gbk
```

### **12. List All Gene Names in the `.gbk` File**
List all gene names using the `/gene=` tag.

```bash
habiba_faried22@cloudshell:~/biocomputing$ grep '/gene=' wildtype.gbk
```

For unique gene names:

```bash
habiba_faried22@cloudshell:~/biocomputing$ grep -o '/gene="[^"]*"' wildtype.gbk | sed 's|/gene="||; s|"||' | sort -u
```

### **13. Clear Your Terminal Space and Print All Commands Used Today**
Clear the terminal and view the history of commands.

```bash
habiba_faried22@cloudshell:~/biocomputing$ clear
habiba_faried22@cloudshell:~/biocomputing$ history
```

### **14. List the Files in Both Folders**
List the files in the `biocomputing` and `habiba` folders.

```bash
habiba_faried22@cloudshell:~/biocomputing$ ls
habiba_faried22@cloudshell:~/habiba$ ls
```

---

## **Project 2: Installing Bioinformatics Software on the Terminal**

### **1. Activate Your Base Conda Environment**
Activate the base environment and check the Conda version.

```bash
habiba_faried22@cloudshell:~$ conda activate base
(base) habiba_faried22@cloudshell:~$ conda --version
```

### **2. Create a Conda Environment Named `funtools`**
Create a new conda environment.

```bash
(base) habiba_faried22@cloudshell:~$ conda create -n funtools
```

### **3. Activate the `funtools` Environment**
Activate the newly created environment.

```bash
(base) habiba_faried22@cloudshell:~$ conda activate funtools
```

### **4. Install Figlet Using Conda or Apt-Get**
Install **Figlet** using `conda` or `apt-get`.

```bash
(funtools) habiba_faried22@cloudshell:~$ conda install -c conda-forge figlet
```

Or using `apt-get`:

```bash
(funtools) habiba_faried22@cloudshell:~$ sudo apt-get update
(funtools) habiba_faried22@cloudshell:~$ sudo apt-get install figlet -y
```

### **5. Run Figlet with Your Name**
Use **Figlet** to display your name.

```bash
(funtools) habiba_faried22@cloudshell:~$ figlet Habiba
```

### **6. Install BWA (Burrows-Wheeler Aligner)**
Install **BWA** through the Bioconda channel.

```bash
(funtools) habiba_faried22@cloudshell:~$ conda install -c bioconda bwa
```

### **7. Install BLAST**
Install **BLAST** using the Bioconda channel.

```bash
(funtools) habiba_faried22@cloudshell:~$ conda install -c bioconda blast
```

### **8. Install Samtools**
Install **Samtools** using the Bioconda channel.

```bash
(funtools) habiba_faried22@cloudshell:~$ conda install -c bioconda samtools
```

### **9. Install Bedtools**
Install **Bedtools** using the Bioconda channel.

```bash
(funtools) habiba_faried22@cloudshell:~$ conda install -c bioconda bedtools
```

### **10. Install Spades**
Install **SPAdes** using the Bioconda channel.

```bash
(funtools) habiba_faried22@cloudshell:~$ conda install -c bioconda spades
```

### **11. Install BCFtools**
Install **BCFtools** using the Bioconda channel.

```bash
(funtools) habiba_faried22@cloudshell:~$ conda install -c bioconda bcftools
```

### **12. Install Fastp**
Install **fastp** using the Bioconda channel.

```bash
(funtools) habiba_faried22@cloudshell:~$ conda install -c bioconda fastp
```

### **13. Install MultiQC**
Install **MultiQC** using the Bioconda channel.

```bash
(funtools) habiba_faried22@cloudshell:~$ conda install -c bioconda multiqc
```

--- 

## **Conclusion**

This repository covers basic shell operations (Bash) and how to install essential bioinformatics tools using **Conda**. The steps provide detailed examples for setting up bioinformatics environments, downloading files, performing file manipulation, and running basic bioinformatics tasks.
