#add all students as a vector
students <- c("Bailey","Ronnie","Ivette", "Areli", "Brittnany","Sarah","Jimbo","Hannah","Silas","Lexie")

#Ranbdomly select students to present
#2 students will present on April 28, 4 on April 30 and 4 on May 5.
# Also, for each presentation we have 4 students with jobs: Scribe, Timekeeper, Moderator, And Assistant.
# We cannot repeat presenters. Each student presents one time. We can repeat jobs, but each student should at least have one job. They should not hold 
# the same job more than once, or more than a job for a single day.
#Set seed for reproducibility
# I will not set seed as it will be random each time you run the code.
#Randomly select students for each presentation day

  April28= sample(students, 2, replace = FALSE)
  April30 = sample(setdiff(students, c(April28)), 4, replace = FALSE)
  May5 = sample(setdiff(students, c(April28, April30)), 4, replace = FALSE)


#Define jobs
jobs <- c("Scribe", "Timekeeper", "Moderator", "Assistant")
#Assign jobs to students for each presentation day
# The students presenting SHOULD NOT have a job on the same day, but they can have a job on a different day.

assign_jobs <- function(presenters) {
  non_presenters <- setdiff(students, presenters)
  assigned_jobs <- sample(jobs, 4, replace = FALSE)
  names(assigned_jobs) <- sample(non_presenters, 4, replace = FALSE)
  return(assigned_jobs)
}

presentation_days <- list(April28, April30, May5)

#Assign jobs for each presentation day
presentation_jobs <- lapply(presentation_days, assign_jobs)
print(presentation_jobs)
