function f=multifun(x,y,normals,idx,varargin)

% parse inputs
n_normals=length(normals);
parser = inputParser;
addRequired(parser,'normals',@isstruct);
addRequired(parser,'idx',@isnumeric);
addParameter(parser,'priors',ones(1,n_normals)/n_normals, @(x) isnumeric(x) && all(x > 0) && all(x < 1));
addParameter(parser,'vals',eye(n_normals), @(x) isnumeric(x) && ismatrix(x));
parse(parser,normals,idx,varargin{:});

vals=parser.Results.vals;
priors=parser.Results.priors;

other_idxs=[1:idx-1, idx+1:n_normals];

xv=[x;y];

% flist=cell(length(other_idxs),1);
% for i=1:length(flist)
%     flist{i}=quad2fun(opt_reg_quad([normals(idx).mu, normals(idx).v],[normals(other_idxs(i)).mu, normals(other_idxs(i)).v],...
%         'prior_1',priors(idx)/(priors(idx)+priors(other_idxs(i))),'vals',vals([idx other_idxs(i)],[idx other_idxs(i)])));
% end

for i=1:n_normals
    maha(i,:)=dot(xv-normals(i).mu,normals(i).v\(xv-normals(i).mu),1);
end

w=exp(-50*maha);
f=sum((1:n_normals)'.*w)./sum(w);

%[~,f]=min(maha);

%f_each=@(x,y) min(cellfun(@(f) f(x,y),flist));
%f=@(x,y) arrayfun(f_each,x,y);
%f=@(x,y) f12(x,y).*f13(x,y);