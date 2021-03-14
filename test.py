from ArtinWedderburn import *
svd_threshold = 0.00001


print("irreps for S2")
ArtinWedderburn(symmetric_group_algebra(2),svd_threshold)
print("\n\n\n")

print("irreps for S3")
ArtinWedderburn(symmetric_group_algebra(3),svd_threshold)
print("\n\n\n")


print("irreps for S4")
ArtinWedderburn(symmetric_group_algebra(4),svd_threshold)
