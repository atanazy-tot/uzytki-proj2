set NURSES ordered;
set SHIFTS ordered;
set DAYS ordered;

param N > 0 integer;
param J > 0 integer;
param K > 0 integer;

param demand {DAYS, SHIFTS};
param vacation {DAYS, NURSES};
param pref_comp {NURSES, NURSES};
param unpref_comp {NURSES, NURSES};
param workhours {NURSES};
param pref_shifts {NURSES, DAYS, SHIFTS};
param unpref_shifts {NURSES, DAYS, SHIFTS};

var x {NURSES, SHIFTS, DAYS} binary;
var z {NURSES, NURSES, SHIFTS, DAYS} binary;

maximize happiness: (sum{i1 in NURSES, i2 in NURSES, j in SHIFTS, k in DAYS} (pref_comp[i1, i2] - unpref_comp[i1, i2])*z[i1, i2, j, k]) + (sum{i in NURSES, j in SHIFTS, k in DAYS} (pref_shifts[i, k, j] - unpref_shifts[i, k, j])*x[i, j, k]);

subject to employer_demand{j in SHIFTS, k in DAYS}: sum{i in NURSES} (1 - vacation[k, i])*x[i, j, k] >= demand[k, j];
subject to max_one_shift_a_day{i in NURSES, k in DAYS}: sum{j in SHIFTS} (1 - vacation[k, i])*x[i, j, k] <= 1;
subject to never_morning_after_night{i in NURSES, k in DAYS: ord(k) < K}: (1 - vacation[k, i])*x[i, member(J, SHIFTS), k] + (1 - vacation[next(k), i])*x[i, member(1, SHIFTS), next(k)] <= 1;
subject to min_working_hours{i in NURSES}: sum{j in SHIFTS, k in DAYS} x[i, j, k] >= workhours[i]*J/24;
subject to max_six_nights_a_week{i in NURSES, k in DAYS: ord(k) mod 7 = 1}: sum{l in ord(k)..(ord(k)+6)} x[i, member(J, SHIFTS), member(l, DAYS)] <= 6; 
subject to justice{i in NURSES}: -2 <= (sum{j in SHIFTS, k in DAYS: ord(k) mod 7 = 0 or ord(k) mod 7 = 6} x[i, j, k]) - (1/N)*(sum{tilde_i in NURSES, j in SHIFTS, k in DAYS: ord(k) mod 7 = 0 or ord(k) mod 7 = 6} x[tilde_i, j, k]) <= 2;
subject to z_conditions1{i1 in NURSES, i2 in NURSES, j in SHIFTS, k in DAYS}: z[i1, i2, j, k] <= x[i1, j, k];
subject to z_conditions2{i1 in NURSES, i2 in NURSES, j in SHIFTS, k in DAYS}: z[i1, i2, j, k] <= x[i2, j, k];
subject to z_conditions3{i1 in NURSES, i2 in NURSES, j in SHIFTS, k in DAYS}: z[i1, i2, j, k] >= x[i1, j, k] + x[i2, j, k] - 1;
subject to enforce_vacation{i in NURSES, j in SHIFTS, k in DAYS}: vacation[k, i]*x[i, j, k] <= 0;