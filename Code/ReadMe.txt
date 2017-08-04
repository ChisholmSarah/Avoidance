This code tests whether dyads avoid each other or seek each others
proximity. This test is based on permutations of time blocks and is
described in detail in the paper Chisholm et al. 2018 (S. Chisholm, A.B.
Stein, N.R. Jordan, T.M. Hubel, J. Shawe-Taylor, T. Fern, J.W. McNutt, A.M.
Wilson, S. Hailes. Parsimonious Test of Dynamic Interaction. Methods in
Ecology and Evolution, 2018).

To run this test use AvoidanceAssociationFun(). An example of how to run
the function is in ExampleScript.m. The detailed description of the input
and output is below.

INPUT:    Id1     - Name of one of the individuals in single quotation
                    marks (for file naming purposes). For example 'John'.
          Id2     - Name of the other individual (for file naming
                    purposes). 
          Data1   - a 4xn1 cell array, where n1 is the number of
                    observations for individual 1. This matrix should
                    contain the location data for the first individual.
                    The first column should be the date ('dd/mm/yyyy' -
                    the single quotation marks are important), the second
                    column should be the time ('HH:MM:SS' - the single
                    quotation marks are important), the third and fourth
                    column are the x and y location in meters from a
                    reference point (which should be the same point for
                    both individuals).
          Data2   - a 4xn2 matrix, where n2 is the number of observations
                    for individual 2. The matrix should have the same
                    format as Data1.
          Dist    - intervals to be examined (Do not include 0 and Inf,
                    the interavals will be 0-Dist(1,1);
                    Dist(1,1)-Dist(1,2); ...; Dist(1,end)-Inf; For
                    example input [50, 100, 150, 200] will test 0-50m,
                    50m-100m, 100m-150m, 150m-200m and 200m-Inf.
          varargin: Use name value pairs, for example if you want to
                  set SigLevel to 0.01, include "'SigLevel', 0.01"
                  after including Id1, Id2, Data1, Data2 and Dist.
              SigLevel    - Significance value to be used for testing.
                            0.05 is the default value. The test adjusts
                            for multiple testing, so this does not need
                            to be taken into acount when choosing the
                            significance level.
              BlockStart  - Time of day the block will start in hours
                            since midnight (for nocturnal animals this
                            could be midday, as individuals are least
                            active during that time). By default this
                            will be midday.
              perm        - number of permutations to be performed. 
                            Default is 10,000.
OUTPUT:   DistpLess   - A matrix with 3 columns. The first column
                        gives the distances (Dist) tested. The
                        second column gives the p-values related to
                        the hypothesis that the two individuals are
                        less often with the tested distances of each
                        other. The third column gives only the
                        p-values for the significant results and NaN
                        for the not significant results.
          DistpMore   - Similarly to DistpLess, but this matrix
                        gives the p-values related to the hypothesis
                        that the individuals are more often within
                        the tested distances of each other.
          p-vlaue plots              - Illustration of the p-values
                                       for the distances (Dist)
                                       tested, related to the
                                       hypothesis that the
                                       individuals were less often
                                       and more often within the
                                       distances of each other than
                                       expected by chance.
          Time series plots          - Time series of the distance
                                       between the 2 individuals.
          Significant location plots - Locations when the two
                                       individuals were within a
                                       significant distance of each
                                       other. The convex hull of the
                                       locations of the two
                                       individuals show a rough
                                       estimate of the individuals'
                                       territories.

For an EXAMPLE run:
Download both folders AvoidanceAssociationFunction and Example and keep 
them in the same folder. Run ExampleScript in the folder Example for
examples.
