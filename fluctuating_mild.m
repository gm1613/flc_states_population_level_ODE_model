function T = fluctuating_mild(t)

day_time = mod(t,1);
t_shift = -13;
time_array = ((0:24))/24;
temperature_array = [5,4,3,3,4,4,4,4,4,6,7,8,8,9,8,7,5,5,5,5,4,4,3,3];

time_shifted_temperature_array = circshift(temperature_array,t_shift);
full_temperature_array = [time_shifted_temperature_array, time_shifted_temperature_array(1)];
T = interp1(time_array,full_temperature_array,day_time);




end