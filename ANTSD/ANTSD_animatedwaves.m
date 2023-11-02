function ANTSD_animatedwaves(rotations)

%% Variables
theta = linspace(0,abs(rotations)*2*pi,50);
r = 1; % radius of circle
X = cos(theta);
if (rotations > 0)
    Y = sin(theta);
else
    Y = -sin(theta);
end

%% Plot

for i = 1:(numel(theta)-1)

    % Ready graphs for this iteration
    pause(0.1);
    cla

    % Plot 1st graph
    subplot(2,2,1)
    % Unit cirlce
    axis equal
    hold on;
    [a,b] = cylinder(r,100);
    plot(a(1,:),b(1,:));
    axis([-2 2 -2 2]) % Set axis
    % x line
    line([0 X(i)],[0 0]);
    % y line
    line([X(i) X(i)],[0 Y(i)]);
    % Hypotenuse line
    line([0 X(i)],[0 Y(i)]);
    % Vertical line
    line([X(i) X(i)],[Y(i) -2],'Color','r');
    % Horizontal line
    line([X(i) 2],[Y(i) Y(i)],'Color','r');
    % Red dot
    plot(X(i),Y(i),'o','MarkerEdgeColor','r',...
        'MarkerFaceColor','r',...
        'MarkerSize',5)

    % Plot 2nd graph
    subplot(2,2,2)
    % sine wave
    if (rotations > 0)
        ezplot( @(theta) sin(theta),[0 theta(i+1)]);
    else
        ezplot( @(theta) -sin(theta),[0 theta(i+1)]);
    end
    axis([0 abs(rotations)*2*pi -2 2]) % Set axis
    % Horizontal line
    if (rotations > 0)    
        line([0 theta(i)],[sin(theta(i)) sin(theta(i))],'Color','r');
    else
        line([0 theta(i)],[-sin(theta(i)) -sin(theta(i))],'Color','r');
    end

    % Plot 3rd graph
    subplot(2,2,3)
    % cos wave
    ezplot( @(theta) cos(theta),[0 theta(i+1)]);
    axis([0 abs(rotations)*2*pi -2 2]) % Set axis
    % Horizontal line
    line([0 theta(i)],[cos(theta(i)) cos(theta(i))],'Color','r');
    view([90 -90])

    % Plot 4th graph
    subplot(2,2,4)
    % sine and cos wave
    if (rotations > 0)
        plot(theta,sin(theta),theta,cos(theta),'Color','b');
    else
        plot(theta,-sin(theta),theta,cos(theta),'Color','b');
    end
    axis([0 abs(rotations)*2*pi -2 2]) % Set axis
    % Horizontal line1
    if (rotations > 0)
        line([0 theta(i)],[sin(theta(i)) sin(theta(i))],'Color','r');
    else
        line([0 theta(i)],[-sin(theta(i)) -sin(theta(i))],'Color','r');
    end
    % Horizontal line2
    line([0 theta(i)],[cos(theta(i)) cos(theta(i))],'Color','r');
end

end