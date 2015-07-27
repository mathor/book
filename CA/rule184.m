clear
rand('seed',0)
N = 100;
ncar = 25;
lane = zeros(1,N);
x = randperm(N,ncar);
lane(x) = 1;
LANE = [];
for i = 1:100
   forward = lane([2:end 1]);
   I = find(lane&~forward);
   lane(I) = 0;
   I(I==N)=0;
   lane(I+1) = 1;
   LANE = [LANE;lane];
end

imagesc(LANE)