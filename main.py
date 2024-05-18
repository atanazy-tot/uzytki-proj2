"""
Main file
"""

from src.input_parser import generate_dat, employer_demands
from src.check_output_sense import check, timetable

input_path = 'data/test_instance.csv'
output_path = 'neos_output_2_variant1.txt'
N, J, K = generate_dat(input_path, 'model/dat_files/first_data.dat')
time, empl_demands = employer_demands(input_path, J, K)

Variables = check(output_path)
timetable(Variables, N, J, K, 'results/archiv/timetable_neos_output_2_variant1.csv', empl_demands)