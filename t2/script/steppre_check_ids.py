def read_ids(file_path):
    with open(file_path, 'r') as f:
        ids = {line.strip() for line in f if line.strip()}
    return ids


file1 = "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/ca2_ids.txt"
file2 = "/storage/caishangLab/fangxiunan/ASAP/test1223_proteincluster_foldseek/ca2_allids_prokaryotes.txt"


ids1 = read_ids(file1)
ids2 = read_ids(file2)


total1 = len(ids1)
total2 = len(ids2)


only_in_file1 = ids1 - ids2
only_in_file2 = ids2 - ids1
common = ids1 & ids2


print("Basic comparison information:")
print("File 1: " + file1)
print("File 2: " + file2)
print("Number of IDs in ca1_ids.txt: " + str(total1))
print("Number of IDs in ca1_allids_prokaryotes.txt: " + str(total2))
print("Number of IDs only in ca1_ids.txt: " + str(len(only_in_file1)))
print("Number of IDs only in ca1_allids_prokaryotes.txt: " + str(len(only_in_file2)))
print("Number of common IDs: " + str(len(common)))