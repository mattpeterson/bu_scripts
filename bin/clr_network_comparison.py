#!/usr/bin/env python
"""
Looks at CLR scores, as calculated from the overexpression microarrays, and
compares them between genes found to have regulations by ChIP and genes which
are not to 
"""
from database import DatabaseConnection
from matplotlib import pyplot as plt
from scipy.stats.distributions import hypergeom

def get_counts(scores, interactions):
    """
    Returns the scores for non-interactions and interactions.  Note that the
    number of interactions for a given gene != number of matches in CLR, since
    there can be doubles of genes.
    """
    int_scores = []
    no_int_scores = []

    # Split CLR scores into those corresponding to interacting genes
    # and those which do not. Need to look for the reverse of the interaction
    # as well, since the CLR score is symmetric.
    for interaction, score in scores.items():
        reverse_interaction = (interaction[1], interaction[0])
        if interaction in interactions:
            int_scores.append(score)
        else:
            if reverse_interaction in interactions:
                # Putting a zero in here, not sure how cool that is.
                no_int_scores.append(0.0)
            else:
                no_int_scores.append(score)

    return int_scores, no_int_scores

def get_clr_probability(scores, int_score, no_int_scores):
    """Calculates a hypergeometric pdf to see how well CLR explains the
    network
    """
    # Get number of CLR "hits" (this threshold should be decided by FDR, 
    # rather than a straight up cutoff)
    int_hits = len(filter(lambda y: y > 1, int_score))
    no_int_hits =len(filter(lambda y: y > 1, no_int_scores))

    # Counts for statistics
    total_hits = int_hits + no_int_hits
    total_genes = len(int_score) + len(no_int_scores)
    num_pulled = len(int_score)

    return 1 - hypergeom.cdf(int_hits, total_genes, total_hits, num_pulled)

def main():
    """Main function"""
    connection = DatabaseConnection()

    # Retrieve interactions from database
    interactions = connection.get_chip_interactions()
    # from_genes = { i[0] for i in interactions }
    from_genes = set([ i[0] for i in interactions ])
    # Retrieve CLR scores from database
    int_scores = []
    no_int_scores = []
    for g in from_genes:
        scores = connection.get_clr_scores(g)
        int_score, no_score = get_counts(scores, interactions)
        int_scores += int_score
        no_int_scores += no_score

        print("%s: %d %f" % (g, len(int_score), 
                            get_clr_probability(scores, int_score, no_score)))

    # Now, plot the distributions
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.set_xlabel('CLR Z-score')
    ax.set_ylabel('Count')
    ax.set_title('Interacting pairs')

    n, bins, patches = ax.hist(int_scores, 25)
    ax = fig.add_subplot(212)
    ax.hist(no_int_scores, bins)
    ax.set_xlabel('CLR Z-score')
    ax.set_ylabel('Count')
    ax.set_title('Non-interacting pairs')

    fig.savefig('CLR_scores.png', dpi=150)

    plt.show()

if __name__ == "__main__":
    main()
