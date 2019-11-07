import pandas as pd
import json
import re

#Read list of pfams
def read_pfam(filename):
    df = pd.read_csv(filename)
    ase_count = 0
    #Use "ase" to distinguish enzymes
    for index, row in df.iterrows():
        if "ase" in row["Pfam Name"]:
            ase_count += 1

    #Number of enzymes and percentage
    print("Ase count =", ase_count)
    print("Percent enzyme =", ase_count/len(df))

    return df

#Classify enzyme words within EC classes
def classify_ECs(filename):
    #EC class lists
    EC_names = [[], [], [], [], [], [], []]
    with open(filename, "r") as f:
        json_f = json.load(f)
        #print(json_f["children"][0]["children"][0]["children"][0]["children"][0]["name"])

        #For loops within for loops
        #1st level - x.
        for i in range(len(json_f["children"])):
            #2nd level - x.x
            for j in range(len(json_f["children"][i]["children"])):
                #3rd level - x.x.x
                for k in range(len(json_f["children"][i]["children"][j]["children"])):
                    #4th level - x.x.x.x
                    for l in range(len(json_f["children"][i]["children"][j]["children"][k]["children"])):
                        #json_f["children"][i]["children"][j]["children"][k]["children"][l]["name"] will find each EC number and Name
                        str = json_f["children"][i]["children"][j]["children"][k]["children"][l]["name"]
                        line = str.split()
                        for word in line:
                            if "ase" in word:
                                #last occurance of a dash
                                ld = word.rfind("-")
                                #append the enzyme name post-dash with no non-alphanumeric characters
                                EC_names[i].append(re.sub("[\W_]+", "", word[ld+1:]))

        #remove repeats and non-ase words
        for i in range(len(EC_names)):
            for ec in EC_names[i]:
                if "ase" not in ec:
                    EC_names[i].remove(ec)
            EC_names[i] = set(EC_names[i])

        return EC_names

def analyze_ECs(EC_names):
    #print distribution of ECs over each class
    # for i in range(len(EC_names)):
    #     print("EC", i+1, "\n")
    #     print(EC_names[i])
    #     print("\n")

    # #Goal = find where overlaps are
    # for i in range(len(EC_names)):
    #     for j in range(len(EC_names)):
    #         if i != j:
    #             print("Intersection of EC", i+1, "and EC", j+1, "is: \n")
    #             print([value for value in EC_names[i] if value in EC_names[j]])
    #             print("\n")

    #Goal - remove overlapping ECs from all but one
    EC_keep = [["monooxygenase", "dehydrogenase", "reductase", "oxygenase", "oxidase", "transhydrogenase", "oxidoreductase"], ["demethylase", "phosphorylase", "cyclotransferase", "glucosyltransferase", "adenylytransferase", "phosphoribosyltransferase", "carboxyltransferase"], ["deaminase", "demethylase", "dehalogenase", "deacetylase", "dextranase", "ribonuclease", "ATPase", "dihydrolase", "diphosphatase", "phosphatase"], ["lyase", "cyclase", "decarboxylase", "hydralase"], ["isomerase", "epimerse"], ["carboxylase", "synthetase", "chelatase", "cobaltochelatase"]]

    #remove all overlapping ECs
    for i in range(len(EC_names)):
        for j in range(len(EC_keep)):
            if i != j:
                EC_names[i] = [value for value in EC_names[i] if value not in EC_keep[j]]

    #return the EC names
    return EC_names

def link_ECs_to_pfams(EC_names, pfam_df):
    for index, row in pfam_df.iterrows():
        #print(index)
        #print(row["Pfam Name"])
        for word in row["Pfam Name"].split():
            if word in EC_names[0]:
                pfam_df.at[index, "EC Class"] = 1
                #row["EC Class"] = 1
                #print(row)
            elif word in EC_names[1]:
                pfam_df.at[index, "EC Class"] = 2
            elif word in EC_names[2]:
                pfam_df.at[index, "EC Class"] = 3
            elif word in EC_names[3]:
                pfam_df.at[index, "EC Class"] = 4
            elif word in EC_names[4]:
                pfam_df.at[index, "EC Class"] = 5
            elif word in EC_names[5]:
                pfam_df.at[index, "EC Class"] = 6
            elif word in EC_names[6]:
                pfam_df.at[index, "EC Class"] = 7
            else:
                pfam_df.at[index, "EC Class"] = 0

    return pfam_df

def EC_pfam_dist(taxonID, pfam_df):
    ECs = []
    pfam_ECs = []
    #EC distribution
    with open("ec_dicts/" + taxonID + "_ecdict.json", "r") as f:
        json_f = json.load(f)
        for e in json_f["enzymes"]:
            ECs.append(int(e[:1]))
        print("EC Distribution:", ECs)

    #pfam EC distribution
    with open("pfam_dicts/" + taxonID + "_pfamdict.json", "r") as f:
        json_f = json.load(f)
        for pfam in json_f["pfams"]:
            pfam_ECs.append(int(pfam_df.loc[pfam_df["Pfam Name"] == pfam]["EC Class"]))
        print("Pfam EC Distribution:", pfam_ECs)



def main():
    #List of all pfams
    pfam_df = read_pfam("pfamlist3582_withEC_07112019.csv")

    #Finding words associated with EC classes
    EC_names = classify_ECs("EC.json")
    EC_names = analyze_ECs(EC_names)

    #Dictionary of pfams - linking pfames to EC numbers
    pfam_df = link_ECs_to_pfams(EC_names, pfam_df)

    #Actual analyses - distribution of enzymes
    EC_pfam_dist("2001200003", pfam_df)


main()
