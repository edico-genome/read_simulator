% These two files appear to contain the same information, but neither appears to be an input template.
% DistMatrix comprises a 4-dimensional matrix:
% 500 x 42 x 4 x 4
% where 500 is the number of cycles = 2 * 250
% 42 is the number of base qualities
% and 4x4 is the reference base x sequenced base transition matrix
%
% QTransMatrix is a 42 x 42 x 498 matrix.
% This apparently gives the probability distribution of a base quality given the
% previous cycle's base quality.  What doesn't make sense to me is that the total
% count fluctuates with each cycle (Eric). This could be due to 'N' bases in the 
% SEQ field not contributing to QTransMatrix counts. The 'N' bases always seem to
% have base_qual=2 ('#') in SRR82646-379_shuffle16k_2.20M.bam.

filename = '/mnt/archive/jcr/reads/NA12878/SRR82646-379.40M/pirs_profile/baseCall_matrix_250bp.count.matrix';
%filename = '/mnt/archive/jcr/reads/NA12878/SRR82646-379.40M/pirs_profile/baseCall_matrix_250bp.ratio.matrix';

fp = fopen(filename, 'r');
done = 0; ctr = 0;
while ~done
    line = fgets(fp);
    if isequal(line,-1), done = 1; 
    else,
       ctr = ctr+1;
       lines{ctr} = line;
    end
end
fclose(fp);


num_qvals = 42;
num_cycles_per_read = 250;

% Do this a different way to keep access to the base at start of line.
for j = 1:4*2*num_cycles_per_read,
    vecs{j} = str2num(['[',lines{11+j}(3:end),']']);
end

mat = cat(1, vecs{:});

vecs = cell(num_qvals,1);
tot_prev_bq_count_per_cycle = zeros(num_qvals,2*(num_cycles_per_read-1));
mat3 = zeros(num_qvals,num_qvals,2*(num_cycles_per_read-1));
% Loop over cycles, getting BQ transition matrix per cycle
for k = 1:2*(num_cycles_per_read-1),
    for j = 1:num_qvals,
        vecs{j} = str2num(['[',lines{2015+j+(k-1)*num_qvals},']']);
    end
    mat2 = cat(1, vecs{:});    
    tot_prev_bq_count_per_cycle(:,k) = mat2(:,end);
    mat2 = mat2(:, 3:end-1);   % Current base qual distribution per preceding base qual val
    mat3(:,:,k) = mat2;
end

% Write header
outdir = '/mnt/archive/jcr/reads/NA12878/SRR82646-379.40M/pirs_profile/synth_pirs_profile_2x500';
outfile = sprintf('%s/baseCall_matrix_2x500bp.count.matrix', outdir);
fid = fopen(outfile,'w');
fprintf(fid,'#Generate @ %s by cooper\n', datestr(now,'HH:MM:SS,yyyy-mm-dd'));
fprintf(fid,'#Input File: [../SRR82646-379_shuffle16k_2.20M.bam]\n');
fprintf(fid,'#Total mapped Reads: 37649156 , mapped Bases 8154645494\n');
fprintf(fid,'#Total statistical Bases: 3633506653 , Reads: 14534757 of ReadLength 250\n');
fprintf(fid,'#Dimensions: Ref_base_number 4, Cycle_number 500, Seq_base_number 4, Quality_number 42\n');
fprintf(fid,'#Mismatch_base: 26345489, Mismatch_rate: 0.725070614037486 %%\n');
fprintf(fid,'#QB_Bases: 267316454, QB_Mismatches: 18191886 (bases with quality <= 2)\n');
fprintf(fid,'#Reference Base Ratio in reads: A 29.831 %%;   C 20.083 %%;   G 19.746 %%;   T 30.34 %%;\n\n');

% Write [DistMatrix] section
fprintf(fid,'[DistMatrix]\n');
fprintf(fid,'#Ref\tCycle');
dna = 'ACGT';
for ii = 1:4,
   for bq = 0:41,
       fprintf(fid,'\t%c-%d',dna(ii),bq);
   end
end
fprintf(fid,'\tRowSum\n');
mi = 0;
for ii = 1:4,
   dbase = dna(ii);
   outcycle = 0;
   for k = 1:2*num_cycles_per_read,
      % Write original line w/ appropriate outcycle
      outcycle = outcycle + 1;
      mi = mi + 1;
      fprintf(fid,'%c\t%d', dbase,outcycle);
      fprintf(fid,'\t%d', mat(mi,:)); fprintf(fid,'\n');
      % Write original line again with next outcycle
      outcycle = outcycle + 1;
      fprintf(fid,'%c\t%d', dbase,outcycle);
      fprintf(fid,'\t%d', mat(mi,:)); fprintf(fid,'\n');
   end
end
fprintf(fid,'<<END\n\n');

% Write [QTransMatrix] section
outcycle = 1;
fprintf(fid,'[QTransMatrix]\n');
fprintf(fid,'#Cycle\tpre1_Q'); fprintf(fid,'\t%d',[0:41]); fprintf(fid,'\tRowSum\n');
for k = 1:2*(num_cycles_per_read-1),
    % write current BQ matrix with new cycle number
    outcycle = outcycle + 1;
    matstr = '';
    for ii = 1:num_qvals,
        prev_qual = ii-1;
        matstr = [matstr sprintf('%d\t%d', outcycle, prev_qual)];
        matstr = [matstr, sprintf('\t%d',mat3(ii,:,k)), sprintf('\t%d\n',tot_prev_bq_count_per_cycle(ii,k))];
    end
    fprintf(fid,'%s',matstr);

    % make diagonal matrix derived from current BQ matrix
    bqtot = sum(mat3(:,:,k),1);  % bqtot is total base count per base qual val for current cycle
    outmat = diag(bqtot);
    outcycle = outcycle + 1;
    matstr = '';
    for ii = 1:num_qvals,
       matstr = [matstr sprintf('%d\t%d', outcycle, ii-1)];
       matstr = [matstr, sprintf('\t%d',outmat(ii,:)), sprintf('\t%d\n',bqtot(ii))];
    end
    fprintf(fid,'%s',matstr);   % write diagonal matrix

    % we may need to write diag matrix again if we're at end of read
    if rem(k,num_cycles_per_read-1) == 0,
        outcycle = outcycle + 1;
        matstr = '';
        for ii = 1:num_qvals,
           matstr = [matstr sprintf('%d\t%d', outcycle, ii-1)];
           matstr = [matstr, sprintf('\t%d',outmat(ii,:)), sprintf('\t%d\n',bqtot(ii))];
        end
        fprintf(fid,'%s',matstr);

        % pirs does not write a BQ matrix for the 1st base of a read, so we skip a cycle
        outcycle = outcycle + 1;
    end
end
fprintf(fid,'<<END\n');
fclose(fid);


if 0,
figure; hold on;
plot(sum(mat3(:,:,1),1), 'b+-');
plot(sum(mat3(:,:,1),2), 'rx-');
plot(sum(mat3(:,:,2),1), 'gd--');
plot(sum(mat3(:,:,2),2), 'mo--');
set(gca,'yscale','log');
end

% We see that sum(mat3(:,:,1),1) is *almost* identical to sum(mat3(:,:,2),2).
% Ignoring the "almost" for now, this suggests that the row represents the base
% quality of the previous cycle, and the column represents the base quality
% of the current cycle.  This seems consistent with the figure above.  We see
% that a previous base quality of 42 can produce a current base quality of 2, but
% not vice versa.  (if the base quality were monotonically decreasing, then
% the matrix would be lower triangular; instead it's merely skewed in that direction).
%figure; imagesc(log10(mat3(:, :, 30))); colorbar

%figure; plot(sum(mat3(:,:,29),1) - sum(mat3(:,:,30),2).', '-o');

% we also see that the total number of elements changes per cycle, and not
% monotonically.
%for j = 1:5, 
%    fprintf('cycle %d, count = %d\n', j, sum(sum(mat3(:,:,j))));
%end


