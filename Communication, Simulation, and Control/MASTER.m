clear s; clear; close all; clc

% Define the serial ports and baud rates
serialPort1 = 'COM6';  % Tachometer data
serialPort2 = 'COM5';  % Voltage, current, and power data
serialPort3 = 'COM4';  % Sensor data
baudRate = 9600;

% Create the serial port objects
s1 = serialport(serialPort1, baudRate);  % Tachometer sensor
s2 = serialport(serialPort2, baudRate);  % Voltage, current, power sensor
s3 = serialport(serialPort3, baudRate);  % Proximity sensor

% Configure terminators
configureTerminator(s1, "CR");
configureTerminator(s2, "CR");
configureTerminator(s3, "CR");

% Initialize variables for storing data
tachometerData = [];
voltageData = [];
currentData = [];
powerData = [];
proximityData = [];
time = [];

% Set up figure for real-time plotting
figure;

% Tachometer data plot
subplot(5,1,1); % Tachometer plot
tachLine = animatedline;
xlabel('Time (s)');
ylabel('Tachometer Data');
title('Real-Time Tachometer Data');
grid on;

% Voltage, Current, Power data plots
subplot(5,1,2); % Voltage plot
voltageLine = animatedline;
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Real-Time Voltage Data');

subplot(5,1,3); % Current plot
currentLine = animatedline;
xlabel('Time (s)');
ylabel('Current (A)');
title('Real-Time Current Data');

subplot(5,1,4); % Power plot
powerLine = animatedline;
xlabel('Time (s)');
ylabel('Power (W)');
title('Real-Time Power Data');

% Proximity sensor data plot
subplot(5,1,5); % Proximity data plot
proximityLine = animatedline;
xlabel('Time (s)');
ylabel('Proximity (mm)');
title('Real-Time Proximity Data');

writeline(s3, "A2");  % Set metric units to mm
writeline(s3, "Z0");  % Define zero point
writeline(s3, "S1000");  % Define sample rate

% Start timer
tic;

while true
    %% Read Proximity Sensor Data
    proximityReading = readline(s3);
    numericProximityData = str2double(proximityReading);
    
    %% Read Tachometer Data
    writeline(s1, '@D0');  % Send command to request tachometer datas
    scrap1 = readline(s1);
    scrap2 = readline(s1);
    tachometerReading = readline(s1);
    numericTachData = str2double(tachometerReading);
    
    %% Read Voltage, Current, Power Data
    writeline(s2, 'GETD');  % Send command to request voltage and current data
    sensorData2 = readline(s2);
    scrap1 = readline(s2);
    voltage = str2double(extractBetween(sensorData2, 1, 4)) / 100;  % Extract and convert voltage
    current = str2double(extractBetween(sensorData2, 5, 8)) / 100;  % Extract and convert current
    power = voltage * current;  % Calculate power

    %% Store the data
    currentTime = toc;  % Get the current time
    
    tachometerData(end+1) = numericTachData;
    voltageData(end+1) = voltage;
    currentData(end+1) = current;
    powerData(end+1) = power;
    proximityData(end+1) = numericProximityData;
    time(end+1) = currentTime;

    %% Update the real-time plots
    addpoints(tachLine, currentTime, numericTachData);  % Update Tachometer plot
    addpoints(voltageLine, currentTime, voltage);  % Update Voltage plot
    addpoints(currentLine, currentTime, current);  % Update Current plot
    addpoints(powerLine, currentTime, power);  % Update Power plot
    addpoints(proximityLine, currentTime, numericProximityData);  % Update Proximity plot
    
    drawnow;
    %pause(0.00112);  % Pause for 0.00125 second
end
