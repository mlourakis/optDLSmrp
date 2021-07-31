% wrapper around Nakano-based MRP solver
% See M. Lourakis and G. Terzakis, "A Globally Optimal Method for the PnP Problem with MRP Rotation Parameterization", ICPR 2020
%
% Arguments: 3xN 3D points, 2xN corresponding calibrated image points
%            Optional arg chk_singular specifies whether the solver should check for a singular matrix and return without a solution in that case
%
% Returns R, t and optionally the corresponding reprojection error

% Manolis Lourakis, Feb. 2020

function [R t minrepr] = optDLS_mrp(p3, p2, chk_singular)

  if(nargin<=2), chk_singular=0; end

  npts=size(p3,2);
	
  % prepare M, see suppl. at http://www.bmva.org/bmvc/2015/papers/paper078/

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %{
	% A=[X1*Y1 ... Xn*Yn];
	% B=[X1 ... Xn];
	A=zeros(2*npts,9);
	B=zeros(2*npts,3);
	for i=1:npts
		X=[0  -1 p2(2,i);  1 0 -p2(1,i)]; % 2x3
		Y=[p3(:,i)' zeros(1,6);  zeros(1,3) p3(:,i)' zeros(1,3); zeros(1,6) p3(:,i)']; % 3x9
		A(2*i-1:2*i,:)=X*Y;
		B(2*i-1:2*i,:)=X;
	end

	%M=A'*A - A'*B*inv(B'*B)*B'*A; % (A-5)
	tMat=inv(B'*B)*B'*A; % (A-3)
	M=A'*A - A'*B*tMat;
  %}
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % more efficient calculation of M; note that the formula for A'*A in the suppl. is wrong!
  Z=zeros(3,3); Zu=zeros(3,3); Zv=zeros(3,3); Zuv=zeros(3,3);
  sump2=zeros(2,1); sump2sq=0;
  up3=zeros(3,1); vp3=zeros(3,1); uvp3=zeros(3,1); sump3=zeros(3,1);
  for i=1:npts
    p2_=p2(:,i);
    p3_=p3(:,i);

    p2_sq=p2_'*p2_;

    % A'*A
    p3p3t=p3_*p3_';
    Z=Z+p3p3t;
    Zu=Zu-p2_(1).*p3p3t;
    Zv=Zv-p2_(2).*p3p3t;
    %Zuv=Zuv + p2_(1)^2.*p3p3t + p2_(2)^2.*p3p3t;
    Zuv=Zuv + p2_sq.*p3p3t;

    % B'*B
    sump2=sump2 + p2_;
    sump2sq=sump2sq + p2_sq;

    % B'*A
    sump3=sump3 + p3_;
    up3=up3 - p2_(1).*p3_;
    vp3=vp3 - p2_(2).*p3_;
    %uvp3=uvp3 + p2_(1)^2.*p3_ + p2_(2)^2.*p3_;
    uvp3=uvp3 + p2_sq.*p3_;
  end
  ze3=zeros(3,3);
  AtA=[Z ze3 Zu; ze3 Z Zv; Zu Zv Zuv];
  BtA=[sump3' 0 0 0 up3';  0 0 0 sump3' vp3';  up3' vp3' uvp3'];
  BtB=[npts 0 -sump2(1); 0 npts -sump2(2); -sump2' sump2sq];

  tMat=inv(BtB)*BtA; % (A-3)
  M=AtA - BtA'*tMat;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % solve
  [b c d]=solver_optDLS_mrp(M, chk_singular);
  nsol=size(b,2);
  if(nsol==0)
    R=[];
    t=[];
    if(nargout>2), minrepr=[]; end
    %warning('optDLS_mrp: no solutions found, singular rotation?');
    return;
  end

  %{
  for i=1:nsol
    psi=[b(i); c(i); d(i)]; % MRP vec
    spsi=-psi./(psi'*psi); % shadow MRP
    for j=i+1:nsol
      if(norm([b(j); c(j); d(j)]-spsi)<1E-04)
        fprintf('solutions %d - %d are equivalent [total %d]\n', i, j, nsol);
      end
    end
  end
  %}

  minerr=Inf;
  minmrp=zeros(3,1);
  minR=zeros(3,3);
  mint=zeros(3,1);
  %cost=zeros(nsol,1);
  for i=1:nsol
    psi=[b(i); c(i); d(i)]; % MRP vec
    magsq=psi'*psi;
    if(magsq>1.0), continue; end % outside unit sphere, will check shadow that is within

    q=[(1-magsq); 2*psi]./(1+magsq); % quat, already normalized
    %R=quat2mat(q); % R
    R=[
      q(1)^2+q(2)^2-q(3)^2-q(4)^2   2*(q(2)*q(3)-q(1)*q(4))       2*(q(2)*q(4)+q(1)*q(3));
      2*(q(2)*q(3)+q(1)*q(4))       q(1)^2-q(2)^2+q(3)^2-q(4)^2   2*(q(3)*q(4)-q(1)*q(2));
      2*(q(2)*q(4)-q(1)*q(3))       2*(q(3)*q(4)+q(1)*q(2))       q(1)^2-q(2)^2-q(3)^2+q(4)^2 ];
    r=reshape(R',[9,1]);
    t=-tMat*r; % (A-3)
    %cost(i)=r'*M*r; % (3)

    xp3=R*p3 + repmat(t, 1,npts);  % [R t]*[p3; ones(1, npts)];
    if(~isempty(find(xp3(3,:) < 0))), continue; end % points must be in front of the camera

    proj=xp3(1:2,:)./repmat(xp3(3,:), 2,1);
    e=proj-p2;
    err=sum(sum(e.^2, 1));
    if(err<minerr)
      minerr=err;
      minmrp=psi;
      minR=R;
      mint=t;
    end
  end
  %fprintf('optDLS_mrp: minimum reprojection error %f for MRP vector [%f %f %f]\n', minerr, minmrp(1), minmrp(2), minmrp(3));
	
  R=minR;
  t=mint;

  if(nargout>2), minrepr=minerr; end
end
