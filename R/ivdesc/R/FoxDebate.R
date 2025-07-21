#' The effects of watching a Fox debate on Proposition 209
#'
#' The data set (n=507) contains findings from the experiment described in Albertson and Lawrence (2009) 
#' in which a representative sample of survey respondents in Orange County, California, were randomly
#' assigned to receive encouragement to view a Fox debate on affirmative action, which would take
#' place on the eve of the 1996 presidential election. Shortly after the election, these respondents were
#' reinterviewed. The postelection questionnaire asked respondents whether they viewed the debate,
#' whether they supported a California proposition (209) to eliminate affirmative action (\code{support}), 
#' and how informed they felt about the proposition (\code{infopro}). The dataset can be used to reproduce 
#' Table 2 in Aronow and Carnegie (2013). Note that mean imputation was used to handle missing data so 
#' non-integer values are imputed. \code{support} and \code{infopro} are excepted and include missing values.
#' 
#' This dataset data documentation has been copied from the archived R package \emph{cicsw}. 
#'
#' @format A data frame with 507 observations on the following 12 variables:
#' \describe{
#'   \item{partyid}{An 11 point scale from "strong Republican" to "strong Democrat".}
#'   \item{pnintst}{Respondent interest in politics and national affairs. Coded 1 = "very interested", 2 = "somewhat interested", 3 = "only slightly interested", 4 = "not interested at all".}
#'   \item{watchnat}{Frequency of national television news consumption. Coded 1 = "never", 2 = "less than once a month", 3 = "once a month", 4 = "several times a month", 5 = "once a week", 6 = "several times a week", 7 = "every day".}
#'   \item{educad}{Education level of respondent. Coded 1 = "eighth grade or less", 2 = "beyond eighth grade, not high school", 3 = "ged", 4 = "high school", 5 = "less than one year vocational school", 6 = "one to two year vocational school", 7 = "two years or more vocational school", 8 = "less than two years of college", 9 = "two or more years of college", 10 = "finished a two-year college program", 11 = "finished a four-year college program", 12 = "master degree or equivalent", 13 = "ph.d., m.d., or other advance degree".}
#'   \item{readnews}{How often respondent reads political news. Coded 1 = "never", 2 = "less than once a month", 3 = "once a month", 4 = "several times a month", 5 = "once a week", 6 = "several times a week", 7 = "every day".}
#'   \item{gender}{Respondent gender. Coded 1 for female and 0 for male.}
#'   \item{income}{Family income from all sources. Coded 1 = "under $10,000", 2 = "between $10,000 and $20,000", 3 = "between $20,000 and $30,000", 4 = "between $30,000 and $40,000", 5 = "between $40,000 and $50,000", 6 = "between $50,000 and $60,000", 7 = "between $60,000 and $70,000", 8 = "between $70,000 and $80,000", 9 = "between $80,000 and $90,000", 10 = "between $90,000 and $100,000", 11 = "$100,000 or more".}
#'   \item{white}{Binary indicator coded 1 if subject is white and 0 otherwise.}
#'   \item{support}{Support for Proposition 209. Coded 1 if subject voted against or opposed and 0 if subject voted for or favored}
#'   \item{infopro}{Information on Proposition 209. Coded from 1 to 4, with 4 meaning respondents had a great deal of information about Proposition 209 prior to the election, and 1 meaning respondents reported no information about the proposition before the election.}
#'   \item{watchpro}{Binary indicator coded 1 if subject watched the Fox Debate about affirmative action and 0 otherwise. This is the outcome ("treatment") of interest.}
#'   \item{conditn}{Binary indicator coded 1 if subject was (randomly) prompted to watch the Fox Debate about affirmative action. This is the encouragement (instrumental) variable.}
#' }
#'
#' @references 
#' Bethany Albertson and Adria Lawrence. (2009). After the credits roll: The long-term effects of educational television on public knowledge and attitudes. \emph{American Politics Research}. 37(2): 275-300.
#' 
#' Peter M. Aronow and Allison Carnegie. (2013). Beyond LATE: Estimation of the average treatment effect with an instrumental variable. \emph{Political Analysis}. 21.4 (2013): 492-506.
#' 
#' Peter M. Aronow and Allison Carnegie. (2013). Replication data for: Beyond LATE: Estimation of the average treatment effect with an instrumental variable. \emph{Dataverse Network.} http://hdl.handle.net/1902.1/21729 (accessed May 14, 2015).
#' 
#' @keywords datasets
#' @docType data
#' @name FoxDebate
"FoxDebate"