context("Multispecies trait models")

test_that("msm defaults produces a glmerMod", {
  msm_glmer <- msm("present", "plot", "logit_rock", "species", "ln_sla", eucs)
  expect_equivalent(class(msm_glmer), "glmerMod")
})
