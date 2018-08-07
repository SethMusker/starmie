# test-clumpp.R
context("checking CLUMPP and CLUMPAK algorithms")

clummp_in <- exampleStructure("clumpp")
# retrieve Q_list
Q_list <- lapply(clummp_in, getQ)
# dummy matrix for IO checking
Q_k1 <- matrix(1L, nrow = nrow(Q_list[[1]]), ncol = 1)

test_that("CLUMPP I/O checks",
          {
            expect_error(clumpp(Q_list, method = "dopey"),
                         "Not a valid CLUMPP method, please use on of: 'greedy', 'greedyLargeK' or 'stephens'")
            expect_error(clumpp(list(Q_list, a = 1:100 )),
                         "cluster runs must be a list of Q matrices")
            expect_error(clumpp(Q_list, iter = -1),
                         "number of iterations must be a positive integer")
            expect_error(clumpp(list(Q_list[[1]], Q_k1)),
                         "size of all matrices in Q_list must be equal")
            expect_identical(clumpp(list(Q_list[[1]])), list(Q_list[[1]])) # checking returns for K=1
            expect_identical(clumpp(list(a = Q_k1, b = Q_k1)), list(a = Q_k1, b = Q_k1)) # checking retruns for multiple K = 1

          })


# test algorithm with K = 6 runs
k6_file <- system.file("extdata/microsat_testfiles/", "locprior_K6.out_f",
                       package = "starmie")
k6_msat <- loadStructure(k6_file)

k6_file_run2 <- system.file("extdata/microsat_testfiles/",
                            "run2_locprior_K6.out_f",
                            package = "starmie")
k6_run2 <- loadStructure(k6_file_run2)
k6_all <- structList(k6_msat, k6_run2)
Q_list <- lapply(k6_all, getQ)

c_greedy <- clumpp(Q_list)
c_largeKgreedy <- clumpp(Q_list, method = "greedyLargeK")

test_that("CLUMMP algorithm checks",
          {
            expect_identical(c_greedy, c_largeKgreedy)
          })
