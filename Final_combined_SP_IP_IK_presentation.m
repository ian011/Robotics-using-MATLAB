% %% Initializing the Recorder Object properties
% % Read a speech signal, and get the speech samples (y) and sampling frequency (fs)
% % recorder = audiorecorder(Fs,nBits,NumChannels);
% clearvars;
% clc;
% 
% 
% Fs = 16000 ; 
% nBits = 8 ; 
% numChannels = 1 ;
% 
% recObj = audiorecorder(Fs,nBits,numChannels);
% 
% %% Recording an audio signal
% disp('Start speaking.')
% recordblocking(recObj, 5);
% disp('End of Recording.');
% 
% %PLaying recorded file
% play(recObj);
% 
% %% Get recorded audio signal in numeric array
% %returns recorded audio data associated with audiorecorder object recorder 
% % in a double array y in the range(-1 to 1)
% Ys = getaudiodata(recObj);
% 
% plot(Ys)
% title('Audio Signal');
% 
% %% Write an Audio File
% % % writes a matrix of audio data, y, with sample rate Fs to a file
% filename = 'speech0.wav';
% audiowrite(filename,Ys,Fs);
% clear y Fs
% 
% %% Reading and Playing back audio
% [Ys,Fs] = audioread(filename);
% sound(Ys,Fs);

%% Create speechClient
% Involves calling the functions from the speech2text folder.
% setting up a scpeechClient object with the speech API and its properties

speechObject = speechClient('Microsoft','language','en-US');

%% Reading in a wav file
[Y, Fs] = audioread('speech20.wav');
% Convert the stereo sound to the mono sound by loading only one channel
% data into the result vector
Ys = Y(:,1);

sound(Ys, Fs)


%% Perform Speech-to-Text Transcription
tableOut = speech2text(speechObject,Ys,Fs,'HTTPTimeOut',25);

 %% Text Extraction from Speech
text = tableOut.Transcript(1) %Returning table data through variable name

preString = regexprep(text, '\W', ''); %  Matches a word not alphanumeric, or underscore
requiredString = regexprep(preString, '_', ''); % Containing one or more underscores.

transcription = split(requiredString, 'and');

expression  = "capture"; % To check if a capture or normal move
matchStr = regexp(transcription(2),expression,'match');
transcription2 = extractBefore(transcription(2), 'capture');

if (matchStr ~= expression)
    [numRows,numCols] = size(transcription);

    assert((numRows == 2 && numCols == 1),'Error! Invalid Statement.')

    condition = ((strlength(transcription(1)) == 2) && (strlength(transcription(2)) == 2));
    if (condition)
        startPos = transcription(1);
        endPos = transcription(2);
    else
        assert(condition,'Error! Command format should be: [Start] and [End]')
    end

else
    [numRows,numCols] = size(transcription);
    [numRows2, numCols2] = size(transcription2);
    
    assert((numRows == 2 && numCols == 1) && (numRows2 == numCols2 == 1),'Invalid! Please input valid Statement.')
    
    condition2 = ((strlength(transcription(1)) == 2) && (strlength(transcription2(1)) == 2));
    if (condition2)
        startPos(1) = transcription2(1);
        endPos(1) = "OUT";
        startPos(2) = transcription(1);
        endPos(2) = transcription2(1);
        
        sequence = [startPos(1) endPos(1) startPos(2) endPos(2)];
    else
        assert(condition,'Error! Command format should be: [Start] and [End] Capture')
    end  
  
end %End of Capture/Normal movement if conditional

%% Text mapping onto Board Coordinates using Voice and board geometry
% boardLength = 20; %In cm
% boardWidth = 20;
% 
% gridDimension = boardLength/8;
% midGrid = gridDimension/2;
% 
% midpointX = 0;
% 
% newMap = containers.Map;
% 
% position = zeros(8); %Not really necessary since storing data in containers
% count = "A";
% for i = 1:8
%     midpointY = 0;
% %     midpointX = midpointX + (midGrid + ((i*2) - 1));
%     midpointX = (midGrid * ((i*2) - 1));
%     coordX = char(i + 64); %Using ASCII Characters;
%     
%     for j = 1:8
%         coordY = num2str(j);
%         coordinates = strcat(coordX, coordY);
%        
% %         midpointY = midpointY + (midGrid + ((j*2) - 1));
%         midpointY = (midGrid * ((j*2) - 1));
%         
%         position = [midpointX midpointY];
%         
%         newMap(coordinates) = position;
%     end
% end
%    

%% Board Grid mapping onto Robot Workspace using image processing baze
% % startPos = 'A1';
% % endPos = 'B2';
% 
% initialCoord = newMap(startPos);
% finalCoord = newMap(endPos);
% 
% combinedCoord = [initialCoord; finalCoord];
% 
% Y_coord = zeros(1, length(combinedCoord));
% X_coord = zeros(1, length(combinedCoord));
% 
% for i = 1: length(combinedCoord)
%     Y_coord(1,i) =  mapfun(combinedCoord(i,2), 0, 6, 11, 17);
%     X_coord(1,i) = mapfun(combinedCoord(i,1), 0, boardLength, 9, -9);
% end
%% 
% Image pre-processing
% 
% 
% reads in an RGB image
J = checkerboard(80, 4, 4) > 0.5;
% I = imread("chess_2.jpg");
[r c]=size(J);
BW = J;


%%


% converts an RGB image to grayscale
% Ig = rgb2gray(I);

% converting to a binary image
%BW = im2bw(Ig, 0.8);
% removing small objects less than 5000 pixels
BW = bwareaopen(BW, 5000);

%montage({I, BW})
%%


% creating a complement of the BW image
BW_comp = imcomplement(BW);

%montage({BW, BW_comp})
%% 
% Detecting the white checkers

% finding the connected objects in a binary image
CC = bwconncomp(BW,4);
L = labelmatrix(CC);
% convert label matrix to RGB
RGB = label2rgb(L);
%imshow(RGB)
%%

% finding connected objects in the complement image
CC_comp = bwconncomp(BW_comp,4);
L_comp = labelmatrix(CC_comp);
% convert label matrix to RGB
RGB_comp = label2rgb(L_comp);

%montage({RGB,RGB_comp})
%% 
% measure properties of the detected objects such as centroid

properties = regionprops(L, 'BoundingBox', 'Area', 'centroid');

% sorting based on the area of the objects detected
[x, idx] = sort([properties.Area], 'descend');
props = properties(idx);
%%

% figure;imshow(RGB)
% hold on;
% for i = 1 : length(props)
% thisBB = props(i).BoundingBox;
% rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
% 'EdgeColor','r','LineWidth',2 )
% plot(props(i).Centroid(:,1), props(i).Centroid(:,2), 'r*')
% text(props(i).Centroid(:,1), props(i).Centroid(:,2), ['(' num2str(props(i).Centroid(:,1)) ',' num2str(props(i).Centroid(:,2)) ')'])
% end
% impixelinfo
%%
properties_comp = regionprops(L_comp, 'BoundingBox', 'Area', 'centroid');

% sorting based on the area of the objects detected
[x, idx] = sort([properties_comp.Area], 'descend');
props_comp = properties_comp(idx);
%%
% 
% figure;imshow(RGB_comp)
% hold on;
% for i = 1 : length(props_comp)
% thisBB = props_comp(i).BoundingBox;
% rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
% 'EdgeColor','r','LineWidth',2 )
% plot(props_comp(i).Centroid(:,1), props_comp(i).Centroid(:,2), 'r*')
% text(props_comp(i).Centroid(:,1), props_comp(i).Centroid(:,2), ['(' num2str(props_comp(i).Centroid(:,1)) ',' num2str(props_comp(i).Centroid(:,2)) ')'])
% end
% impixelinfo
%%
% while true
% imshow(RGB)
% hold on;
% for i = 1 : length(props)
% thisBB = props(i).BoundingBox;
% rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
% 'EdgeColor','r','LineWidth',2 )
% plot(props(i).Centroid(:,1), props(i).Centroid(:,2), 'r*')
% text(props(i).Centroid(:,1), props(i).Centroid(:,2), ['(' num2str(props(i).Centroid(:,1)) ',' num2str(props(i).Centroid(:,2)) ')'])
% end
% hold on;
% for i = 1 : length(props_comp)
% thisBB = props_comp(i).BoundingBox;
% rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
% 'EdgeColor','r','LineWidth',2 )
% plot(props_comp(i).Centroid(:,1), props_comp(i).Centroid(:,2), 'r*')
% text(props_comp(i).Centroid(:,1), props_comp(i).Centroid(:,2), ['(' num2str(props_comp(i).Centroid(:,1)) ',' num2str(props_comp(i).Centroid(:,2)) ')'])
% end
% impixelinfo
%     
%     break
% end
%% 
% Combining everything into one.
% 
% 
comb = {};
if max(props_comp(1).Centroid>props(1).Centroid) == 0
for i = 1:length(props)
   comb = [comb, props_comp(i).Centroid];
   comb = [comb, props(i).Centroid];
end
else
   for i = 1:length(props)
   comb = [comb, props(i).Centroid];
   comb = [comb, props_comp(i).Centroid];
   end
    
end
comb = comb';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comb_ar = cell2mat(comb);
comb_ar(:,2) = repmat(unique(comb_ar(:,2)), length(comb_ar)/length(unique(comb_ar(:,2))), 1);

for i = 1:64
    checker_bd{i} = [comb_ar(i,1), comb_ar(i,2)];
end
checker_bd = checker_bd';

%% 
% 
% Plotting the chessboard with coordinates
% imshow(RGB)
% hold on;
% for i = 1 : length(comb_ar)
% plot(comb_ar(i,1), comb_ar(i,2), 'r*')
% text(comb_ar(i,1), comb_ar(i,2), ['(' num2str(comb_ar(i,1)) ',' num2str(comb_ar(i,2)) ')'])
% end
% impixelinfo
%% 
% Mapping pixel values to board coordinates

for i = 1:length(comb_ar)
    real_coord{i} =  [mapfun(comb_ar(i,1), 0, r, 9, -9),mapfun(comb_ar(i,2), 0, 192, 11, 17)];
end
real_coord = real_coord';


%  creating key value pairs for board letter and number
key = [];
for char ='A':'H' 
    
    for num = 1:8
        
        key = [key, {append(char, num2str(num))}];
        
    end
    
    
    
end
key = key';

keySet_r = key;
valueSet_r = real_coord;
mapping_r = containers.Map(keySet_r,valueSet_r);

mapping_r('OUT') = [13 11.25]; 



combinedCoord =[];
if (matchStr ~= expression)
    initialCoord = mapping_r(startPos);
    finalCoord = mapping_r(endPos);
    
    combinedCoord = [initialCoord; finalCoord];
    
else
    initialCoord = zeros(length(startPos));
    finalCoord = zeros(length(endPos));
    for i = 1:2
        initialCoord(i,:) = mapping_r(startPos(i));
        finalCoord(i,:) = mapping_r(endPos(i));
        
        combinedCoord = [combinedCoord; initialCoord(i,:); finalCoord(i,:)];
   end
end

Y_coord = zeros(1, length(combinedCoord));
X_coord = zeros(1, length(combinedCoord));

for i = 1: length(combinedCoord)
    Y_coord(1,i) =  combinedCoord(i,2);
    X_coord(1,i) = combinedCoord(i,1);
end

%% Manipulator handing over (FK AND IK)
 if ~isempty(instrfind)
          fclose(instrfind);
          delete(instrfind);
 end

arduino=serial('COM9','Baudrate',9600); %created an arduino object
fopen(arduino); %connects the serial port object to device
pause(2);

% % 
 startup_rvc;
%VOICE INPUT
% Z_in = 6;
% Y_in = Y_coord;
% X_in = X_coord;

%TRAJECTORY ARRAY!!
% 
Z_in = [3.5];
X_in = [-3.375  1.125  5.625  -5.625];
Y_in = [ 12.0417  14.125  16.2083   16.2083];



points = size(Y_in,2); %Getting the number of columns from array
%disp(points);

syms theta1 theta2 theta3;

% % arduino=serial('COM8','Baudrate',9600); %created an arduino object
% 
% fopen(arduino); %connects the serial port object to device
% pause(2);
theta1 = 90 ; %initial position angle configurations 
theta2 = 90;
theta3 = 90;

for init = 1:points
    x = X_in(init);
    y = Y_in(init);
    z = Z_in;
    
    assert((-9 < x) && (x < 9) ,'Error! point outside workspace.');  % defining dextrous workspace
    assert((9 < y) && (y < 18),'Error! point outside workspace.')  ; % defining dextrous workspace
    assert((-9 < z) && (z < 9),'Error! point outside workspace.') ; % defining dextrous workspace
    
   fk_speech = round([x, y, z]);
    
    % Inverse Kinematics Start
    
%%mode=input("Kindly enter the desired mode:   " );
%%if mode==2

%     xyz=input("Kindly enter the desired xyz co-ordinates:   ","s");
%     xyzint=str2num(xyz)
%     x=xyzint(1)
%     y=xyzint(2)
%     z=xyzint(3)

initial_theta1 = theta1;  %initial position angle configurations for trajectory tracking
initial_theta2 = theta2;
initial_theta3 = theta3;

q1_pos = [initial_theta1 initial_theta2 initial_theta3];
    
       
    if x>=0
        theta1=round(atand(x/y)+90);
    else
        thet1=round(atand(y/x));
        theta1=abs(thet1);
    end




    thet3=acosd(((x^2+y^2+(z-5)^2-8^2-12^2)/(2*8*12))); % 8 and 12 are the 
    % dimensions of the upper arm and bottom arm. 
    
    theta3=round((180-thet3)+0);
    
    theta2=round(atand((sqrt(x^2+y^2))\(5.5-z))+atand((8+12*cosd(thet3))\(12*sind(thet3))));
    if theta2 > 0
        theta2 = 90 - abs(theta2) + 90;
    end
    
    anglin=[theta1 theta2 theta3];
    
    assert((theta1 == real(theta1)) && (theta2 == real(theta2)) && (theta3 == real(theta3)),'Error! one of the angles is complex.')%defining dextrous workspace
 q2_pos = [theta1 theta2 theta3];
 t = [0:0.05:2];
 q = mtraj(@tpoly, q1_pos, q2_pos, t);%using a quintic polymomial trajectory generator
 
 imidd_theta1 = round(q(:,1));
 imidd_theta2 = round(q(:,2));
 imidd_theta3 = round(q(:,3));
 
 Manip_imidd_thetas = [imidd_theta1 imidd_theta2 imidd_theta3];
 
 figure,
 plot(t, q); %desired actual motion
 title("Ideal Motion" + num2str(init));
 
%  thetaexperiment1 = num2str(imidd_theta1)+"A" ;
%  thetaexperiment1 = join(transpose(thetaexperiment1));
%  thetaexperiment1 =  strtrim( thetaexperiment1 );

%  thetaexperiment2 = join(transpose(num2str(imidd_theta2)+"B"));
 

%  fprintf(arduino,thetaexperiment1);
%   
%  pause(2.0);
 
for expt = 1:length(imidd_theta1)
    
thetaexperiment1 = num2str(imidd_theta1(expt)) + "A" ;
 thetaexperiment1 = join(transpose(thetaexperiment1));
 thetaexperiment2 = join(transpose(num2str(imidd_theta2(expt)) + "B"));
 thetaexperiment3 = join(transpose(num2str(imidd_theta3(expt)) + "C"));
 theta_expt_total = thetaexperiment1 + thetaexperiment2 + thetaexperiment3;
 
 fprintf(arduino,theta_expt_total);
 pause(0.05);
 
end
% %Reading Data from Arduino
% for i = 1:41
%     data=string(fscanf(arduino));
%     datas(i,1) = data;
% end


figure,
 plot(t, Manip_imidd_thetas); %desired actual motion
 title("Actual Manipulator  Motion" + num2str(init));
 
% var = fscanf(arduino);
% disp(var);
%  flushinput(arduino);
%  flushoutput(arduino);
%  fwrite(arduino,thetaexperiment2);
%  pause(2.0);
%  
%  fwrite(arduino,thetaexperiment3);
%  pause(2.0);
%  
%  flushinput(arduino);

%  
 
%  f =  fscanf(arduino, '%s');
%  
% 
% 
%  if ~isempty(instrfind)
%           fclose(instrfind);
%           delete(instrfind);
%    end

%   %  fprintf(arduino,theta1);
%     %output=fscanf(arduino);
%     %disp(output);
%     
    theta4 = 180;
    
%     thetatotal=join([theta1 "A" theta2 "B" theta3 "C"])
%     thetatotal = num2str(theta1)+"A"+num2str(theta2)+"B"+num2str(theta3)+"C"+num2str(theta4)+"D";
%     thetatotal
%     fprintf(arduino,thetatotal);
%     
%     pause(3.0);
% %     
% %  f =  fscanf(arduino, '%s');
%  
% 
%     
%      %forward kinematics start
%    %FK FOR FEEDBACK
% syms delta1;
% syms delta2;
% syms delta3;
% 
% 
% 
% % Find numbers in the cell 
% numbers_cell = regexp(f, '\d*', 'match');
% 
% delta_angles = str2double(numbers_cell);
% delta1 = delta_angles(1);
% delta2 = delta_angles(2);
% delta3 = delta_angles(3);
% 
% %delta4 = delta_angles(4);
% fwd_k_input = [delta1 delta2 delta3];
% % INSERT FK CODE HERE
% 
% syms th d alph a;
% %% from sprong book 
% A = trotz(th)*transl(0,0,d)*transl(a,0,0)*trotx(alph);
% 
% l1=-5.5; %offset
% l2=0;
% l3=8.0;
% l4=12.0;
% 
% %Link1
% syms th1 L1
% A1=subs(A,{a,alph,d,th},{l2,pi/2,l1,th1});
% 
% %Link2
% syms th2 L2
% A2=subs(A,{a,alph,d,th},{l3,0,0,th2});
% 
% %Link3
% syms th3 L3
% A3=subs(A,{a,alph,d,th},{l4,0,0,th3});
% 
% At=simplify(A1*A2*A3);
% fk=At(:,4);
% 
% %angles = [71 143 117];
% angles = [delta1 delta2 delta3];
% angles(1,3)=180-(angles(1,3));
% %angles(1,2)=180-(angles(1,3))
% angles_rad=angles*pi/180;
% 
% 
% Ep=subs(fk,{th1,th2,th3},{angles_rad});
% Ep1=vpa(Ep,2);
% Ep1(2,1)=(Ep1(2,1)*-1);
% Ep1(3,1)=(Ep1(3,1)*-1);
% Ep1 = round(Ep1); 
% fwd_conf = double(Ep1');
% % t = [0:0.05:2];
% % mdl_puma560;
% % T1 = SE3(x, y, z) * SE3.Ry(pi/2);
% % T2 = SE3(0.5, -0.3, 0.44) * SE3.Ry(pi/2);
% % Ts = ctraj(T1, T2, length(t));
% 
% % % qc = p560.ikine6s(Ts);
% 
% if isequal(fk_speech,fwd_conf)
% 
%     disp("Arduino Execution Complete!")
% else 
%     disp("Execution Error")
% end
% 
% 
% 
% 
% 
% 
% 
%     
% %     if mod(init,2) == 1 %Finding the modulus
% %         fprintf(arduino,"180D90B90C90A"); %Close Gripper
% %         pause(0.5);
% %         
% %     else
% %         fprintf(arduino,"120D90B90C90A"); %Open Gripper
% %         pause(0.5);
% %     end 
% %     
%     
%     
%   
% % %FK END (FEEDBACK)
end


function output = mapfun(value, fromLow, fromHigh, toLow,toHigh)
%Function to map the angle input by user to a suitable angle rotation angle
narginchk(5,5)
nargoutchk(0,1)
output = (value - fromLow) .* (toHigh - toLow) ./ (fromHigh - fromLow) + toLow;
end