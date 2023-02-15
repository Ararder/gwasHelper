test_that("slurm_header can handle variable args", {
  expect <-  c(
    "#!/bin/bash",
    "#SBATCH --time=24:00:00",
    "#SBATCH --mem=24gb",
    "#SBATCH --output=/test_dir/test/slurm-%j.out",
    "\n"
    )

  expect2 <-  c(
    "#!/bin/bash",
    "#SBATCH --time=24:00:00",
    "#SBATCH --mem=24gb",
    "#SBATCH --output=/test_dir/test/slurm-%j.out",
    "#SBATCH --dependency=afterok:100",
    "\n"
  )


  expect_equal(
    expect,
    slurm_header(time = 24, mem =24, output_dir = "/test_dir/test") %>% as.vector()
  )

  expect_equal(
    expect2,
    slurm_header(time = 24, mem =24, dependency = 100, output_dir = "/test_dir/test") %>% as.vector()
  )
})




