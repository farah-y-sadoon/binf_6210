##_ Fork and Clone ------

# this will create a new project in RStudio in a new window

usethis::create_from_github(
  "https://github.com/evainnocente/assignment2-6210.git",
  # replace the above github url with the green Code > Local > Clone using the web URL
  # you need to update this for each repository you want to contribute to
  destdir = "./", 
  # or whatever directory location you want
  # this destdir can be the same across different repositories if you want all your github projects to be in the same location
  # each repository will be in a different subdirectory in this location
  fork = TRUE
)

# Checking set-up
usethis::git_remotes() # Should see an origin in my github account and upstream (Eva's github account)
usethis::git_sitrep()

# Create new branch
usethis::pr_init(branch = "feature/6210_assignment/farahs_edits") # change formidable w/ more informative branch name!

usethis::git_sitrep() # Making sure the current upstream branch is the one I just created

usethis::pr_push()

# Create pull request in browser, then do the following
usethis::pr_pull() # make sure branch is up-to-date
usethis::pr_finish() # cleaning everything up, deleting local 
usethis::git_sitrep()
