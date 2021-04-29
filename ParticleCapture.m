%% Particle Tracking for Capture

% This code requires two tiff stack files of the exact same data, one
% of these files is the raw tiff stack as captured on the Chronos camera.
% And the other is a processed tiff stack.

% The processing procedure is as follows:

% 1) Read the stack of tiffs into Fiji (ImageJ)
        % a) Make a copy by file->Save As and save it as a tiff and give it
        % an appropriate name, I usually put a P at the end of the file
        % name to identify the file as the processed data.
% 2) Use the oval tool in ImageJ to make a circle right around the edge
%    of the drum so that all particles and the chain are included, but the
%    background is not.  Clear the background by edit->Clear Outside. This
%    will make the background black.
% 3) Apply a gaussian blur to the image.  To do this select
%    Process->Filters->Gaussian Blur.  And run a blur of size 7.
% 4) Threshold the image.  Select Image->Adjust->Threshold and the
%    Thresholding window should show up.  Select the dark background and do
%    not reset range checkboxes, and then select the set button.  Set the
%    low threshold to 5 and the high threshold to 35, as long as lighting
%    is the same and the Chronos aperture is completely open.  This should
%    result in just the chain showing up as red.  Then press apply and only
%    select the black background box.
%    DOUBLE CHECK THIS
% 5) There will probably be a ring of white around the edge of the drum
%    after this, and we just want to see the chain so draw another circle
%    around just the inside of the white rim and clear the outside, the
%    resulting image will be completely black with a white outline of the
%    chain.

%% Step 1: Read in the first Image of Both Files

% The original unprocessed file
o = double(imread('A1O.tif',1));

% The processed file
p = im2uint8(imread('A1P.tif',1));

% How many frames
f = length(imfinfo('A1O.tif'));

% Frame dimensions
w = length(p(:,1));
h = length(p(1,:));


%% Step 2: Find the particles in the first frame of the raw image

% Apply bpass filter, the first entry is the raw image, the second should
% almost always be 1, and the third should be an estimate of the size of
% the particles in pixels.  For more info read the bpass.m file.
b = bpass(o,1,15);

% Next we identify the features that we are interested in using the
% pkfnd.m function.  The first entry is the filtered image 'b', the next
% is the threshold value (particles have pixel values above the threshold
% value), and the third is the diameter of the particles.
pk = pkfnd(b,40,15);

% Then we find the centroid of the particles using the cntrd.m function.
cnt = cntrd(b,pk,15);

% Now make a scatterplot of the particles found to check that it is working
% well. NOTE: the 1152 number corresponds to the dimension of the image.
figure(1)
scatter(cnt(:,1),h-cnt(:,2),'r.')
ylim([0,h])
xlim([0,w])
set(gca,'YTick',[])
set(gca,'XTick',[])

%% Now we will find the convex hull of the chain for the first frame.

% First we make the image binary, it really is already binary because of
% the thresholding in ImageJ, but we want the 255 values to be 1 in order
% to run a black and white convex hull function.

pb = p;

for i=1:w
    for j=1:h
        if p(i,j) == 255
            pb(i,j) = 1;
        end
    end
end

% Then find the convex hull of the chain.
ph = bwconvhull(pb);

% Now check that this works
figure(2)
colormap('gray'), imagesc(ph);
ylim([0,h])
xlim([0,w])
set(gca,'YTick',[])
set(gca,'XTick',[])

%% Determine which partiles are captured

% Round particle positions to integers (pixels), and put the positions into
% a matrix that has a third column of all zeros, this collumn we will use 
% as an identifier for captured particles, captured particles will have 
% a value of 1 in the third column.

rp = [round(cnt(:,1)) round(cnt(:,2)) zeros(length(cnt(:,1)),1)];

% Put a value of 1 in the third column for each captured particle.
for i = 1:length(rp(:,1))
    rp(i,3) = ph(rp(i,2),rp(i,1));
end

% Isolate the captured particles and their positions.

cap = find(rp(:,3)==1);
cappart = [rp(cap,1) rp(cap,2)];

% Plot the convex hull and the captured particles together

figure(3)
colormap('gray'), imagesc(ph);

hold on

scatter(cappart(:,1),cappart(:,2),'g.')
ylim([0,h])
xlim([0,w])
set(gca,'YTick',[])
set(gca,'XTick',[])

hold off

% If all is working then we should see the convex hull and green dots for
% only the particles that are captured and inside the convex hull.  If 
% this is not the case, something went wrong, PANIC.

%% Next, we want to track the particles through the entire dataset

% We will make a scramble matrix where the first 2 columns are positions
% of each particle and the third column is the frame number.

scramble = [cnt(:,1) cnt(:,2) ones(length(cnt(:,1)),1)];

% Now loop through for each frame
for j=2:f % Loop to the size of the stack
    o = double(imread('A1O.tif',j));
    b = bpass(o,1,15);
    pk = pkfnd(b,40,15);
    cnt = cntrd(b,pk,15);
    add2scramble = [cnt(:,1) cnt(:,2) (ones(length(cnt(:,1)),1)+j-1)];
    scramble = [scramble;add2scramble];
end

% Note matlab will not be happy that scramble is not preallocated for size
% matlab can suck it up

% Now tie the particles together with the function track.m

t = track(scramble,10);

% Plot trajectories to check that everything is working
figure(4)
q = find(t(:,4)==1);

plot(t(q,1),h-t(q,2))
hold on
for k=2:max(t(:,4))
    q = find(t(:,4)==k);
    plot(t(q,1),h-t(q,2))
end
hold off
ylim([0,h])
xlim([0,w])
set(gca,'YTick',[])
set(gca,'XTick',[])
title('Particle Trajectories')

%% Now determine capture at each frame

nc = zeros(f,1);
figure(5)
hold on
%phs = zeros(w,h,f);
info = [];
%phs(:,:,1) = ph;
for j=1:f % Loop to the size of the stack
    p = im2uint8(imread('A1P.tif',j));
    pb = p;
    for i=1:w
        for k=1:h
            if p(i,k) == 255
            pb(i,k) = 1;
            end
        end
    end
    %phs(:,:,j) = bwconvhull(pb);
    ph = bwconvhull(pb);
    q = find(t(:,3)==j);
    tt = [t(q,:) zeros(length(q),1)];
    %pos = [round(tt(:,1)) round(tt(:,2)) zeros(length(tt(:,1)),1)];
    pos = round(tt);
    for i = 1:length(pos(:,1))
        pos(i,5) = ph(pos(i,2),pos(i,1));
    end
    cap = find(pos(:,5)==1);
%     nocap = find(pos(:,3)~=1);
%     cappart = [pos(cap,1) pos(cap,2)];
%     nocappart = [pos(nocap,1) pos(nocap,2)];
    nc(j) = length(cap);
%     scatter(cappart(:,1),h-cappart(:,2),5,'r.')
%     scatter(nocappart(:,1),h-nocappart(:,2),5,'g.')
    info = [info;pos];
end
hold off
% The vector 'nc' tells us how many particles are captured at each frame
% The matrix 'info' has a is the rounded trajectories matrix with a fifth
% column that has a value of 1 for captured particles and a value of 0 for
% non-captured particles.

%% Save info

csvwrite('info.csv',info)


%% Manipulating info matrix to make interesting plots

d = find(info(:,5)==1);
evercap = unique(info(d,4));

info = [info zeros(length(info(:,1)),1)];

for m = 1:length(info(:,1))
    for n = 1:length(evercap)
        dummy = evercap(n);
        if info(m,4)==dummy
            info(m,6) = 1;
        end
    end
end
%% Plot

ac = find(info(:,5)==1);
c = find(info(:,5)==0 & info(:,6)==1);
nev = find(info(:,6)==0);

%% Break
figure(6)
hold on
% for i = 1:length(info(:,1))
%     if info(i,5)==1 && info(i,6)==1
%         scatter(info(i,1),h-info(i,2),5,'r.')
%     elseif info(i,5)==0 && info(i,6)==1
%         scatter(info(i,1),h-info(i,2),5,'b.')
%     else
%         scatter(info(i,1),h-info(i,2),5,'y.')
%     end
% end

scatter(info(ac,1),h-info(ac,2),1,'r.')
scatter(info(c,1),h-info(c,2),1,'b.')
scatter(info(nev,1),h-info(nev,2),1,'g.')



ylim([0,h])
xlim([0,w])
set(gca,'YTick',[])
set(gca,'XTick',[])
title('Particle Trajectories')
legend('Captured Particles','Non-Captured Particles','Never Captured Particles')
hold off
%% Plot number of captured particles

figure(7)


tiledlayout(2,1)
va = readtable('Values.csv');
% time = 1:length(va(:,1));
% value = va(:,2);

nexttile
stackedplot(va(:,2))
title('Average Pixel Differences')
xlabel('Time')
%ylabel('Average Pixel Value')

nexttile
scatter(1:length(nc),nc,'r.')
% set(gca,'YTick',[])
% set(gca,'XTick',[])
title('Captured Particle Count')
xlabel('Time')
ylabel('Captured Particles')