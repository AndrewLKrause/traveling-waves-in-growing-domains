
% Plot a kymograph
f = figure;
f.Position(1:2) = f.Position(1:2)*0.2; f.Position(3:4) = f.Position(3:4)*2;
Ts = repmat(ts,1,N);

    if(centred == 1)
        Xs = xs - xs(:,end)/2;
    else
        Xs = xs;
    end

surfir(Xs(:),Ts(:),us(:),1);
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$t$', 'Interpreter', 'latex')
c = colorbar;
%colormap(viridis)
c.TickLabelInterpreter = 'latex';
c.Label.String = '$u$';
c.Label.Interpreter = 'latex';
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
axis tight;
view(0,90)
box on

%hold on
%plot(X,X/2,'r','linewidth',2)

function h=surfir(x,y,z,s,opt)

%% default parameters

if (nargin<4)||isempty(s)                                                   % no shrink factor provided
    s=0.5;                                                                   % default value
end

if nargin<5                                                                 % no options provided
    opt={'FaceColor','interp','edgecolor','none'};                           % default
end


%% Remove duplicate data points

[xy,ind] = unique([x,y],'rows');
z=z(ind);
x=xy(:,1);
y=xy(:,2);


%% triangulate data points

dt = delaunayTriangulation(x, y);                                           % Delaunay triangulation

x=dt.Points(:,1);
y=dt.Points(:,2);

%% find enclosing boundary

k=boundary(x,y,s);                                                          % define boundary enclosing all data points
c=[k(1:end-1),k(2:end)];                                                    % constraints

dt.Constraints=c;                                                           % add constraints
io = dt.isInterior();                                                       % triangles which are in the interior of domain
tri=dt.ConnectivityList;                                                    % read triangles
tri=tri(io,:);                                                              % use only triangles in the interior

%% plot

h=trisurf(tri,x,y,z-2,z,opt{:});                                              % plot triangles and interpolate colors
end