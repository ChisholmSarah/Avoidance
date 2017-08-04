%% Run this sript for examples

%% Example without Association
load NoAssociation_D50_OAD150_IAD74_AT5.mat
cd ../AvoidanceAssociationFunction/

[DistpLessNoAssociation, DistpMoreNoAssociation] = AvoidanceAssociationFun('ID1', 'ID2', NoAssociation1, NoAssociation2, [200, 400, 600, 800], 'SigLevel', 0.05, 'BlockStart', 12, 'perm', 1000);

%% Example with Association
cd ../Example/
load Association_D50_OAD150_IAD74_AT5.mat
cd ../AvoidanceAssociationFunction/

[DistpLessAssociation, DistpMoreAssociation] = AvoidanceAssociationFun('ID1', 'ID2', Association1, Association2, [200, 400, 600, 800], 'SigLevel', 0.05, 'BlockStart', 12, 'perm', 1000);
