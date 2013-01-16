import argparse
import csv
import MotifSearch

parser = argparse.ArgumentParser(description='Analyze data files generated by MotifSearch iterators.')
parser.add_argument('data_file')
parser.add_argument('-c')
parser.add_argument('-i')
args = parser.parse_args()

num_header_rows = 4

# Compare dataset to original motif
if args.c:
    original_motif = MotifSearch.read_plain(args.c)
    with open(args.data_file) as data_file:
        csv_reader = csv.reader(data_file)
        [csv_reader.next() for i in range(num_header_rows)]  # Eliminate header rows
        iterations_conserved_sites = []
        iteration_size = 0
        for row in csv_reader:
            i = int(row[0])
            if i is 0:
                iteration_size += 1
            if row[1] in original_motif:
                try:
                    iterations_conserved_sites[i] += 1
                except:
                    iterations_conserved_sites.append(1)
    print "Iteration   Conservation (%) | Iteration   Conservation (%)"
    for i, iteration in enumerate(iterations_conserved_sites):
        cell = "%9d   %d/%d (%.2f%%)" % (i, iteration, iteration_size, float(iteration) / iteration_size * 100)
        print cell.ljust(29) + "|",
        if i % 2 is 1 and i is not 0:
            print
