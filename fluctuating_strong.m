function T = fluctuating_strong(t)

day_time = mod(t,1);
t_shift = -13;
time_array = ((0:24))/24;
temperature_array = [-1,-1,-1,1,1,2,2,4,6,9,11,11,12,12,11,10,10,7,4,4,2,2,2,0];

time_shifted_temperature_array = circshift(temperature_array,t_shift);

full_temperature_array = [time_shifted_temperature_array, time_shifted_temperature_array(1)];
T = interp1(time_array,full_temperature_array,day_time);




end