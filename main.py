"""
Main file
"""

from input_parser import *
from check_output_sense import *

input_path = 'test_instance.csv'
output_path = 'neos_output_2_variant1.txt'
N, J, K = generate_dat(input_path, 'first_data.dat')
time, empl_demands = employer_demands(input_path, J, K)

Variables = check(output_path)
timetable(Variables, N, J, K, 'timetable_neos_output_2_variant1.csv', empl_demands)